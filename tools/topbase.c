/* ------- file: -------------------------- readtopbase.c ----------- */

/* Reads atomic data compilation from TopBase
 *
 * Han Uitenbroek
 * Last modified: Fri Apr 29 03:50:24 2011 --
 */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "rh.h"
#include "atom.h"
#include "constant.h"
#include "error.h"
#include "inputs.h"

#define NLAMB_DEFAULT   15
#define QCORE_DEFAULT  1.0
#define QWING_DEFAULT 30.0
#define VDW_DEFAULT    1.0


/* --- Function prototypes --                          -------------- */

int  countTokens(char *line, char *separator);
void strupcase(char *string);
void writeModelAtom(Atom *atom, FILE *fpOut);


/* --- Global variables --                             -------------- */

CommandLine commandline;
char messageStr[MAX_LINE_SIZE];


/* ------- begin -------------------------- readTopBase.c ----------- */

int main(int argc, char *argv[])
{
  register int i, j, la;
  
  char    input[MAX_LINE_SIZE], parity[2], config1[10], config2[10],
          state[4], orbital[2], spectrumID[5];
  bool_t  exit_on_EOF, readphoto = FALSE, hunt, lines = FALSE;
  int     Nlevel, Nl, Nread, itb, *islp, *ilv, ln, ll, multiplicity,
          islp_p, ilv_p, i_p, nz, ne, nprint, Nlamb,
          stage, i_l, islp_l, ilv_l, jslp_l, jlv_l, Nrad, ntoken;
  double  gi, *te, e, rl, abundance, weight, Eion, *alpha,
         *E, hERc, Econt, fdum, gf, f, lambda0, C, Ecorr, E0;
  FILE   *fp_energy, *fp_photo, *fp_lines;
  Atom atom;
  AtomicLine *line;
  AtomicContinuum *continuum;

  commandline.quiet = FALSE;
  commandline.logfile = stderr;

  if (argc <= 2) {
    fprintf(stderr, "Usage: %s energyFile Nlevel", argv[0]);
    fprintf(stderr, " [photoFile linesFile]\n");
    exit(0);
  }
  Nlevel = atoi(argv[2]);
  if (argc >= 4) readphoto = TRUE;
  if (argc >= 5) lines = TRUE;
  
  hERc = (HPLANCK / E_RYDBERG) * CLIGHT;
  
  /* --- Open the data files --                        -------------- */
  
  if ((fp_energy = fopen(argv[1], "r")) == NULL) {
    sprintf(messageStr, "Unable to open input file %s", argv[1]);
    Error(ERROR_LEVEL_2, argv[0], messageStr);
  }
  if (readphoto) {
    if ((fp_photo = fopen(argv[3], "r")) == NULL) {
      sprintf(messageStr, "Unable to open input file %s", argv[3]);
      Error(ERROR_LEVEL_2, argv[0], messageStr);
    }
  }
  if (lines) {
    if ((fp_lines = fopen(argv[4], "r")) == NULL) {
      sprintf(messageStr, "Unable to open input file %s", argv[4]);
      Error(ERROR_LEVEL_2, argv[0], messageStr);
    }
  }
  /* --- Read atom ID, spectral designation, abundance,
     atomic weigtht, and Energy of the continuum -- ----------------- */
  
  fprintf(stderr, "Give atom ID (char[2]) > "); 
  Nread = scanf("%2s", atom.ID);
  strupcase(atom.ID);
  fprintf(stderr, "Give ionization stage (int) > ");
  Nread = scanf("%d", &stage);
  fprintf(stderr, "Give abundance relative to Hydrogen > ");
  Nread = scanf("%lf", &abundance);
  fprintf(stderr, "Give atomic weight [amu] > ");
  Nread = scanf("%lf", &weight);
  fprintf(stderr, "Give energy of the continuum [E_RYDBERG] > ");
  Nread = scanf("%lf", &Econt);
  fprintf(stderr, "Give energy of the lowest level [E_RYDBERG] > ");
  Nread = scanf("%lf", &E0);
  
  /* --- Count the actual number of levels in the energy table -- --- */
  
  Nl = 0;
  while (getLine(fp_energy, "=", input, exit_on_EOF=FALSE) != EOF) Nl++;
  atom.Nlevel = MIN(Nlevel+1, Nl);
  islp = (int *) malloc(atom.Nlevel * sizeof(int));
  ilv  = (int *) malloc(atom.Nlevel * sizeof(int));
  rewind(fp_energy);

  atom.Nline = atom.Ncont = atom.Nfixed = atom.Nprd = 0;
  atom.weight = weight;
  atom.abundance = (abundance > 0.0) ? abundance : 1.0E-12;
  
  switch (stage) {
  case 0: strcpy(spectrumID, "I");    break;
  case 1: strcpy(spectrumID, "II");   break;
  case 2: strcpy(spectrumID, "III");  break;
  case 3: strcpy(spectrumID, "IV");   break;
  case 4: strcpy(spectrumID, "V");    break;
  default: ;
  }
  
  atom.E = (double *) malloc(atom.Nlevel * sizeof(double));
  atom.g = (double *) malloc(atom.Nlevel * sizeof(double));
  atom.label = (char **) malloc(atom.Nlevel * sizeof(char *));
  atom.stage = (int *) malloc(atom.Nlevel * sizeof(int));
  
  te = (double *) malloc((atom.Nlevel - 1) * sizeof(double));
  
  for (i = 0;  i < atom.Nlevel-1;  i++) {
    getLine(fp_energy, "=", input, exit_on_EOF=TRUE);
    ntoken = countTokens(input, " ");
    switch (ntoken) {
    case 10:
      sscanf(input, "%d %d %d %d %d %s %lf %lf %lf %lf",
	     &itb, islp+i, ilv+i, &ln, &ll,
	     config1, &gi, te+i, &e, &rl);  break;
    case 11:
      sscanf(input, "%d %d %d %d %d %s %s %lf %lf %lf %lf",
	     &itb, islp+i, ilv+i, &ln, &ll,
	     config1, config2, &gi, te+i, &e, &rl);  break;
    case 12:
      sscanf(input, "%d %d %d %d %d %s %s %s %lf %lf %lf %lf",
	     &itb, islp+i, ilv+i, &ln, &ll, 
	     config1, config2, state, &gi, te+i, &e, &rl);  break;
    default: ;
    }
    atom.g[i] = gi;
    atom.E[i] = (te[i] + E0) * E_RYDBERG;
    atom.stage[i] = stage;
    
    if (islp[i] % 2) strcpy(parity, "O");  else  strcpy(parity, "E");
    multiplicity = islp[i] / 100;
    switch ((islp[i] - multiplicity*100) / 10) {
    case  0: *orbital = 'S';  break;
    case  1: *orbital = 'P';  break;
    case  2: *orbital = 'D';  break;
    case  3: *orbital = 'F';  break;
    case  4: *orbital = 'G';  break;
    case  5: *orbital = 'H';  break;
    case  6: *orbital = 'I';  break;
    case  7: *orbital = 'K';  break;
    case  8: *orbital = 'L';  break;
    case  9: *orbital = 'M';  break;
    case 10: *orbital = 'N';  break;
    case 11: *orbital = 'O';  break;
    case 12: *orbital = 'Q';  break;
    case 13: *orbital = 'R';  break;
    case 14: *orbital = 'T';  break;
    case 15: *orbital = 'U';  break;
    default: ;
    }
    *(orbital+1) = '\0';
    
    atom.label[i] = (char *) calloc(ATOM_LABEL_WIDTH+1, sizeof(char)); 
    switch (ntoken) {
    case 10:
      nprint = sprintf(atom.label[i], "%s %s %s %1d%s%s",
		       atom.ID, spectrumID, config1,
		       multiplicity, orbital, parity);  break;
    case 11:
      nprint = sprintf(atom.label[i], "%s %s %s %s %1d%s%s",
		       atom.ID, spectrumID, config1, config2,
		       multiplicity, orbital, parity);  break;
    case 12:
      nprint = sprintf(atom.label[i], "%s %s %s %s %s %1d%s%s",
		       atom.ID, spectrumID, config1, config2, state,
		       multiplicity, orbital, parity);  break;
    default: ;
    }
    strupcase(atom.label[i]); 
  }
  /* --- Finally, the continuum level --               -------------- */
  
  atom.g[atom.Nlevel-1] = 1.0;
  atom.E[atom.Nlevel-1] = (Econt + E0) * E_RYDBERG;
  atom.label[atom.Nlevel-1] =
    (char *) calloc(ATOM_LABEL_WIDTH+1, sizeof(char));
  strcpy(atom.label[atom.Nlevel-1], "Continuum");
  atom.stage[atom.Nlevel-1] = stage + 1;
  
  /* --- Now the radiative transitions if necessary -- -------------- */
  
  Nrad = (atom.Nlevel - 1) * (atom.Nlevel - 2) / 2 + atom.Nlevel - 1;
  atom.line = (AtomicLine *) malloc(Nrad * sizeof(AtomicLine));
  
  /* --- Find the appropriate bound-bound transitions -- ------------ */
  
  atom.Nline = 0;
  line = atom.line;
  if (lines) {
    C = 2 * PI * (Q_ELECTRON/EPSILON_0) * (Q_ELECTRON/M_ELECTRON) / CLIGHT;

    while (getLine(fp_lines, "=", input, exit_on_EOF=FALSE) != EOF) {
      Nread = sscanf(input, "%d %d %d %d %d %lf %lf",
		     &i_l, &islp_l, &ilv_l, &jslp_l, &jlv_l, &gf, &fdum);

      for (i = 0;  i < atom.Nlevel;  i++) {
	if ((islp[i] == islp_l)  &&  (ilv[i] == ilv_l)) {
	  for (j = 0;  j < atom.Nlevel;  j++) {
	    if ((islp[j] == jslp_l)  &&  (ilv[j] == jlv_l)) {
	      if (gf < 0.0) {
		f = -gf / atom.g[i];
		line->i = i;  line->j = j;
	      } else {
		f = gf / atom.g[j];
		line->i = j;  line->j = i;
	      }

	      lambda0 = (HPLANCK * CLIGHT) /
		(atom.E[line->j] - atom.E[line->i]);
	      line->Aji = C / SQ(lambda0) *
		(atom.g[line->i] / atom.g[line->j]) * f;
	      line->Bji = CUBE(lambda0) / (2.0 * HPLANCK * CLIGHT) * line->Aji;
	      line->Bij = (atom.g[line->j] / atom.g[line->i]) * line->Bji;
	      line->lambda0 = lambda0 / NM_TO_M;
	      line->Nlambda = NLAMB_DEFAULT;

	      line->Voigt = TRUE;
	      line->PRD = FALSE;
	      line->cvdWaals[2] = line->cvdWaals[0] = VDW_DEFAULT;
              line->Grad = line->cStark = 0.0;
	      line->qcore = QCORE_DEFAULT;
	      line->qwing = QWING_DEFAULT;

	      line++;  atom.Nline++;
	    }
	  }
	}
      }
    }
  }
  /* --- Find the bound-free transitions --            -------------- */
  
  atom.Ncont = 0;
  atom.continuum =
    (AtomicContinuum *) malloc((Nrad - atom.Nline)* sizeof(AtomicContinuum));
  continuum = atom.continuum;
  if (readphoto) {
    while (getLine(fp_photo, "=", input, exit_on_EOF=FALSE) != EOF) {

      Nread = sscanf(input, "%d %d %d %d %d %lf %d",
		     &i_p, &nz, &ne, &islp_p, &ilv_p, &Eion, &Nlamb);
      E     = (double *) malloc(Nlamb * sizeof(double));
      alpha = (double *) malloc(Nlamb * sizeof(double));
      for (la = 0;  la < Nlamb;  la++)
	Nread = fscanf(fp_photo, "%lf %lf", E+la, alpha+la);

      for (i = 0;  i < atom.Nlevel;  i++) {
	if ((islp[i] == islp_p)  &&  (ilv[i] == ilv_p)) {
	  continuum->hydrogenic = FALSE;
	  continuum->i = i;
	  continuum->j = atom.Nlevel-1;
	  continuum->Nlambda = Nlamb;
	  
	  continuum->lambda =
	    (double *) malloc(continuum->Nlambda * sizeof(double));
	  continuum->alpha =
	    (double *) malloc(continuum->Nlambda * sizeof(double));

	  /* --- Correct the ionization energies towards the energy
	         values given in the table file --     -------------- */

	  Ecorr = Econt - te[i] + Eion;
	  for (la = 0;  la < continuum->Nlambda;  la++) {
	    continuum->lambda[la] = hERc /
	      ((E[(continuum->Nlambda - 1) - la] + Ecorr) * NM_TO_M);
	    continuum->alpha[la]  = alpha[(continuum->Nlambda - 1) - la] *
	      MEGABARN_TO_M2;
	  }
	  /* --- Determine the wavelength and crossection at the true
                 treshold --                           -------------- */

	  Eion = Econt - te[i];
	  continuum->lambda0 = hERc / (Eion * NM_TO_M);
	  Linear(continuum->Nlambda, E, alpha,
		       1, &Eion, &(continuum->alpha0), hunt=FALSE);
	  continuum->alpha0 *= MEGABARN_TO_M2;

	  continuum++;  atom.Ncont++; 	  
	  break;
	}
      }
      free(E);  free(alpha);
    }
  }
  atom.line =
    (AtomicLine *) realloc(atom.line, atom.Nline * sizeof(AtomicLine));
  atom.continuum = (AtomicContinuum *)
    realloc(atom.continuum, atom.Ncont * sizeof(AtomicContinuum));

  /* --- Write the atomic model to file --             -------------- */

  writeModelAtom(&atom, stdout);
}
/* ------- end ---------------------------- readTopBase.c ----------- */

/* ------- begin -------------------------- strupcase.c ------------- */

void strupcase(char *string)
{
  while (*string) {
    *string = toupper(*string);
    string++;
  }
}
/* ------- end ---------------------------- strupcase.c ------------- */

/* ------- begin -------------------------- countTokens.c ----------- */

int countTokens(char *line, char *separator)
{
  char *tmpStr;
  int   count = 1;

  tmpStr = (char *) malloc((strlen(line) + 1) * sizeof(char));
  strcpy(tmpStr, line);
  strtok(tmpStr, separator);
  while (strtok(NULL, separator)) count++;

  free(tmpStr);
  return count;
}
/* ------- end ---------------------------- countTokens.c ----------- */
