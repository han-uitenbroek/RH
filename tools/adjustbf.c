/* ------- file: -------------------------- adjustbf.c -------------- */

/* Reads atomic data file and resamples boundfree data
   according to script file.
 *
 * Han Uitenbroek
 * Last modified: Thu Mar 11 15:10:15 1999 --
 */

#include <malloc.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "constant.h"
#include "error.h"
#include "inputs.h"

/* --- Function prototypes --                          -------------- */

void lowPassFilter(int N, double *afft, double cutoff);
void rawAtom(struct Atom *atom, char *atomFileName);
void realft(int N, double *data, int sign);
void writeModelAtom(struct Atom *atom, FILE *fpOut);

/* --- Global variables --                             -------------- */


/* ------- begin -------------------------- adjustBF.c -------------- */

void main(int argc, char *argv[])
{
  register int kr, la;

  char    label[ATOM_LABEL_WIDTH + 1], keyword[20], inputLine[MAX_LINE_SIZE];
  bool_t  hunt, exit_on_EOF, present;
  int     Nft, Nlamb, nread;
  double   hERc, *E, *lambda, *alpha, *Eft, *lft, dE, dlamb, fdum,
          lambdaMin, cutoff;
  double *aft;
  FILE   *fp_script;
  struct  Atom atom;
  struct  RadTrans *cont;

  hERc = (HPLANCK / E_RYDBERG) * CLIGHT;

  if (argc <= 2) {
    fprintf(stderr, "Usage: %s atomFile scriptFile \n", argv[0]);
    exit(0);
  }
  /* --- Open the data files --                        -------------- */
  
  rawAtom(&atom, argv[1]);
  if ((fp_script = fopen(argv[2], "r")) == NULL)
    Error(UNABLE_TO_OPEN_INPUTFILE, argv[0], argv[2]);

  /* --- Go through the lines of the script --         -------------- */

  while (getLine(fp_script, "*", inputLine, exit_on_EOF=FALSE) != EOF) {
    nread = sscanf(inputLine, "'%20c' %s %d %f",
		   label, keyword, &Nlamb, &fdum);

    present = FALSE;
    for (kr = atom.Nline;  kr < (atom.Nline + atom.Ncont);  kr++) {
      cont = atom.radtrans + kr;  
      if (strstr(label, atom.label[cont->i])) {
	present = TRUE;
	if (strstr(keyword, "FFT")) {

	  E = (double *) malloc(cont->Nlamb * sizeof(double));
	  for (la = 0;  la < cont->Nlamb;  la++) {
	    E[la] = hERc / (cont->lambda[la] * NM_TO_M);
	    cont->alpha[la] /= MEGABARN_TO_M2;
	  }
	  /* --- Use Fourier transform in energy space to get rid
		 of fine structure --                  -------------- */

	  cutoff = fdum;
	  Nft = 1;
	  while (Nft < cont->Nlamb) Nft <<= 1;

	  Eft = (double *) malloc(Nft * sizeof(double));
	  aft = (double *) malloc(2*Nft * sizeof(double));
	  dE = (E[cont->Nlamb-1] - E[0]) / (Nft - 1);
	  Eft[0] = E[0];
	  for (la = 1;  la < Nft;  la++)  Eft[la] = Eft[la-1] + dE;
	  Linear(cont->Nlamb, E, cont->alpha,
		       Nft, Eft, aft, hunt=TRUE);

	  /* --- Symmetrize to avoid aliasing --       -------------- */

	  for (la = 0;  la < Nft;  la++)
	    aft[Nft+la] = aft[Nft-1 - la];

	  /* --- Transform, filter, and transform back -- ----------- */

	  realft(Nft, aft, 1);
	  lowPassFilter(Nft, aft, cutoff);
	  realft(Nft, aft, -1);

	  lft = Eft;
	  for (la = 0;  la < Nft;  la++) {
	    lft[la] = hERc / (Eft[la] * NM_TO_M);
	    aft[la] = fabs(aft[la]);
	  }
	  /* --- Reinterpolate back on equidistant grid -- ---------- */

	  cont->alpha  =
	    (double *) realloc(cont->alpha, Nlamb*sizeof(double));
	  cont->lambda =
	    (double *) realloc(cont->lambda, Nlamb*sizeof(double));

	  dlamb = (lft[0] - lft[Nft-1]) / (Nlamb - 1);
	  cont->lambda[0] = lft[0];
	  for (la = 1;  la < Nlamb;  la++)
	    cont->lambda[la] = cont->lambda[la-1] - dlamb;

	  Linear(Nft, lft, aft, Nlamb,
		       cont->lambda, cont->alpha, hunt=TRUE);
	  for (la = 0;  la < Nlamb;  la++)
	    cont->alpha[la] *= MEGABARN_TO_M2 / Nft;
	  cont->Nlamb = Nlamb;

	  free(E);  free(Eft);  free(aft);
	} else if (strstr(keyword, "HYDROGENIC")) {

	  /* --- Hydrogenic approximation is valid --   ------------- */

	  lambdaMin = fdum;
	  cont->shape = BF_HYDROGENIC;

	  cont->alpha  =
	    (double *) realloc(cont->alpha, Nlamb*sizeof(double));
	  cont->lambda =
	    (double *) realloc(cont->lambda, Nlamb*sizeof(double));

	  cont->Nlamb = Nlamb;
	  dlamb = (cont->lambda0 - lambdaMin) / (Nlamb - 1);
	  cont->lambda[0] = lambdaMin;
	  for (la = 1;  la < Nlamb;  la++)
	    cont->lambda[la] = cont->lambda[la-1] - dlamb;
	  for (la = 0;  la < cont->Nlamb;  la++) {
	    cont->alpha[la] = cont->alpha0 *
	      CUBE(cont->lambda[Nlamb-1] / cont->lambda[la]);
	  }
	} else if (strstr(keyword, "RESAMPLE")) {

	  /* --- Resample wavelengths --              --------------- */

	  lambdaMin = MAX(fdum, lambda[0]);
	  lambda = (double *) malloc(cont->Nlamb * sizeof(double));
	  alpha  = (double *) malloc(cont->Nlamb * sizeof(double));

	  for (la = 0;  la < cont->Nlamb;  la++) {
	    lambda[la] = cont->lambda[la];
	    alpha[la]  = cont->alpha[la];
	  }
	  cont->alpha  =
	    (double *) realloc(cont->alpha, Nlamb*sizeof(double));
	  cont->lambda =
	    (double *) realloc(cont->lambda, Nlamb*sizeof(double));

	  dlamb = (lambda[cont->Nlamb-1] - lambdaMin) / (Nlamb - 1);
	  cont->lambda[0] = lambdaMin;
	  for (la = 1;  la < Nlamb;  la++)
	    cont->lambda[la] = cont->lambda[la-1] + dlamb;
	  Linear(cont->Nlamb, lambda, alpha,
		       Nlamb, cont->lambda, cont->alpha, hunt=TRUE);
	  
	  cont->Nlamb = Nlamb;
	  free(lambda);  free(alpha);
	}
      }
    }
    if (!present)
      fprintf(stderr, "Could not find matching energy level for label %s\n",
	      label);
  }
  writeModelAtom(&atom, stdout);
}
/* ------- end ---------------------------- adjustBF.c -------------- */
