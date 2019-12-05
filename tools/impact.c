/* ------- file: -------------------------- impact.c ----------------

       Version:       rh1.0, tools
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Tue May  9 13:19:23 2006 --
       Last Update/JB Wed Apr 12 15:42:04 MET DST 2000 
       --------------------------                      ----------RH-- */

/* --- Reads atomic model from raw atomic data file atomFile, computes
       collisional cross sections, and prints them to file outFile.

       Use Seatons impact parameter approximation for optically allowed
       atomic bound-bound transitions
  See: Seaton, M. J. 1962, Proc. Phys. Soc. 79, 1105-1117.

       Van Regemorter's approximation is used for ionic bound-bound and
       for atomic bound-bound transitions for which the oscillator
       strength is not given in the atomic data file.

  See: Landolt-Boernstein, Volume 2, Astronomy and Astrophysics,
        subvolume b, Stars and Star clusters, Springer-Verlag, 98-100).

       In this latter case a predifined f-value F_QUADRUPOLE is used.
       --                                              --------------- */
 
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "constant.h"
#include "error.h"
#include "inputs.h"


#define  F_QUADRUPOLE  0.1
/*#define  NTEMP         7 */
#define  NTEMP         13
#define  R_AVG_FLAG    1
#define  VERBOSE       0

#define  COMMENT_CHAR  "#"


/* --- Function prototypes --                          -------------- */

double Bessel_K0(double x);
double Bessel_K1(double x);
void   crossSection(Atom *atom, int Ntemp, double *temp, FILE* fp_out);
double findbeta1(double righthand, double *zeta1);
double GaussLaguerre(double (*function) (double x));
char **getWords(char *label, char *separator, int *count);
void   impactParam(AtomicLine *line, bool_t flag,
		   int Ntemp, double *temp, double *CE);
double phiImpact(double x);
void   quantumNumbers(Atom *atom);
void   rawAtom(Atom *atom, char *atomFileName);
double square(double x) { return x*x; }


/* --- Global variables --                             -------------- */

int    *n, *l;
double *n_eff, E_Rydberg, f, deltaE, R0, x0, ERkT, PIa0square;

CommandLine commandline;
char messageStr[MAX_LINE_SIZE];


/* ------- begin -------------------------- impact.c ---------------- */

int main(int argc, char *argv[])
{
  int    Ntemp = NTEMP;
  /*  double temp[NTEMP] = {3.0E+3, 5.0E+3, 7.0E+3, 1.0E+4, 2.0E+4,
      3.0E+4, 1.0E+5}; */

  double temp[NTEMP] = {1264.9, 2249.4, 4000.0, 7113.1, 12649.1, 22493.6, 40000.0,
			71131.1, 126491.0, 224936.0, 400000.0, 711312.0, 1264910.0};   
  FILE  *fp_out;
  Atom atom;

  commandline.quiet = FALSE;
  commandline.logfile = stderr;

  if (argc >=3) {
    if ((fp_out = fopen(argv[2], "w")) == NULL) {
      sprintf(messageStr, "Unable to open output file %s", argv[1]);
      Error(ERROR_LEVEL_2, argv[0], messageStr);
    }
  } else if (argc == 2)
    fp_out = stdout;
  else {
    fprintf(stderr, "  Usage: %s atomFile [outFile]\n", argv[0]);
    exit(0);
  }
  rawAtom(&atom, argv[1]);
  E_Rydberg = E_RYDBERG / (1.0 + M_ELECTRON / (atom.weight * AMU));

  quantumNumbers(&atom);	  
  crossSection(&atom, Ntemp, temp, fp_out);
  fclose(fp_out);
}
/* ------- end ---------------------------- impact.c ---------------- */

/* ------- begin -------------------------- quantumNumbers.c -------- */

void quantumNumbers(Atom *atom)
{
  register int i;

  char  multiplet[ATOM_LABEL_WIDTH + 1], *ptr, **words;
  int   ic, count;
  double Z;

  /* --- Determine the principal quantum number n, the orbital quantum
         number l, and the effective quantum number n_eff for all
         levels from the atomic data labels. --        -------------- */

  n = (int *) malloc((atom->Nlevel - 1) * sizeof(int));
  l = (int *) malloc((atom->Nlevel - 1) * sizeof(int));
  n_eff = (double*) malloc((atom->Nlevel - 1) * sizeof(double));

  for (i = 0;  i < atom->Nlevel-1;  i++) {

    /* --- The principal and orbital quantum numbers n and l -- ----- */

    strcpy(multiplet, atom->label[i]);
    ptr = multiplet + (strlen(multiplet) - 1);
    while ((*ptr != 'E')  &&  (*ptr != 'O')  &&  (ptr > multiplet))  ptr--;
    if (ptr > multiplet)
      *(ptr + 1) = '\0';
    else {
      sprintf(messageStr, "Unable to determine whether level is even or odd %s",
	      atom->label[i]);
      Error(ERROR_LEVEL_2, "quantumNumbers", messageStr);
    }

    words = getWords(multiplet, " ", &count);
    n[i]  = (int) (words[count-2][0] - '0');
    switch (words[count-1][1]) {
    case 'S': l[i] = 0;  break;
    case 'P': l[i] = 1;  break;
    case 'D': l[i] = 2;  break;
    case 'F': l[i] = 3;  break;
    case 'G': l[i] = 4;  break;
    case 'H': l[i] = 5;  break;
    case 'I': l[i] = 6;  break;
    default: l[i] = n[i] - 1;
    }

    /* --- The effective quantum number n_eff --       -------------- */

    ic = i + 1;
    while ((atom->stage[ic] < atom->stage[i]+1) && (ic < atom->Nlevel))
      ic++;
    if (atom->stage[ic] == atom->stage[i]) {
      sprintf(messageStr, "Found no overlying continuum for level %s",
	      atom->label[i]);
      Error(ERROR_LEVEL_2, "quantumNumbers", messageStr);
    } else {
      Z = (double) (atom->stage[i] + 1);
      n_eff[i] = Z * sqrt(E_Rydberg / (atom->E[ic] - atom->E[i]));
    }
    free(words);
  }
  /* --- Print quantum numbers to standard error --    -------------- */

  fprintf(stderr, "\n       label              n_eff    n   l\n");
  fprintf(stderr, "-----------------------------------------\n");
  for (i = 0;  i < atom->Nlevel-1;  i++) {
    fprintf(stderr, "'%20s'  %f  %2d  %2d\n",
	    atom->label[i], n_eff[i], n[i], l[i]);
  }
  fprintf(stderr, "\n");
}
/* ------- end ---------------------------- quantunNumbers.c -------- */

/* ------- begin -------------------------- getWords.c -------------- */

char **getWords(char *label, char *separator, int *count)
{
  char **theWords;
  int    length = strlen(label);

  /* --- Get the separate words that constitute the label -- -------- */

  *count = 1;
  theWords = (char **) malloc((length/2 + 1) * sizeof(char *));
  theWords[0] = strtok(label, separator);
  while ((theWords[*count] = strtok(NULL, separator)))
    *count += 1;

  return (char **) realloc(theWords, *count * sizeof(char *));
}
/* ------- end ---------------------------- getWords.c -------------- */

/* ------- begin -------------------------- crossSection.c ---------- */

void crossSection(Atom *atom, int Ntemp, double *temp, FILE *fp_out)
{
  register int i, j, k, kr, kf;

  char    methodStr[40];
  bool_t  validtransition;
  int     ic;
  double  gbar_i;
  double *CE, *CI, *Omega, C, C0, C1, C2_atom, C2_ion, C3, alpha0;
  AtomicLine *line;
  AtomicContinuum *continuum;
  FixedTransition *ft;

  /* --- Some useful constants --                      -------------- */

  C  = 2 * PI * (Q_ELECTRON/EPSILON_0) * (Q_ELECTRON/M_ELECTRON) / CLIGHT;
  C0 = 1.55E+11;                            /* --- [s^-s m K^1/2] --- */
  C1 = ((E_RYDBERG/sqrt(M_ELECTRON)) * PI*SQ(RBOHR)) *
    sqrt(8.0/(PI*KBOLTZMANN));              /* --- [s^-1 m^3 K^1/2] - */
  C2_atom = 2.15E-6;                        /* --- [s^-1 m^3 K^3/2] - */
  C2_ion  = 3.96E-6;                        /* --- [s^-1 m^3 K^3/2] - */
  C3 = sqrt(8 * KBOLTZMANN / (PI * M_ELECTRON));
                                            /* --- [s^-1 m K^-1/2] -- */
  PIa0square = PI * square(RBOHR);          /* --- [m^2] --        -- */

  /* --- Prints the atomic collision rate coefficients for bound-bound
         and bound-free transitions of atom pointed to by atomic
         data structure atom --                        -------------- */

  /* --- Allocate space for rate coefficients --       -------------- */

  CE = CI = Omega = (double *) calloc(Ntemp, sizeof(double));

  /* --- Print used temperatures --                    -------------- */

  fprintf(fp_out, "\n%1s Collisional rate coefficients from impact.c\n\n",
	  COMMENT_CHAR);
  fprintf(fp_out, " TEMP   %2d       ", Ntemp);
  for (k = 0;  k < Ntemp;  k++)
    fprintf(fp_out, "%9.1f%2s", temp[k], (k == Ntemp-1) ? "\n\n" : "  ");

  for (i = 0;  i < atom->Nlevel;  i++) {
    for (j = i + 1;  j < atom->Nlevel;  j++) {

      /* --- Only select transitions within the same stage -- ------- */

      if (atom->stage[j] == atom->stage[i]) {
	deltaE = atom->E[j] - atom->E[i];

        /* --- Use either the oscillator strength from the atomic 
               data file, or if that is not present, a predefined
               value F_QUADRUPOLE --                   -------------- */

        validtransition = FALSE;
	for (kr = 0;  kr < atom->Nline;  kr++) {
	  line = atom->line + kr;
	  if (line->i == i  &&  line->j == j) {
	    validtransition = TRUE;
	    f = line->Aji * (atom->g[j]/atom->g[i]) *
	      square(line->lambda0*NM_TO_M) / C;
	    break;
	  }
	}
	for (kf = 0;  kf < atom->Nfixed;  kf++) {
	  ft = atom->ft + kf;
	  if (ft->type == FIXED_LINE  &&  ft->i == i  &&  ft->j == j) {
	    validtransition = TRUE;
	    f = ft->strength;
	    break;
	  }
	}
	if (!validtransition) f = F_QUADRUPOLE;
	
	if (atom->stage[i] == 0) {

          /* --- Transition within neutral stage --    -------------- */

	  fprintf(fp_out, " CE     %2d %2d    ", j, i);
	  if (validtransition) {
	    impactParam(atom->line + kr, R_AVG_FLAG, Ntemp, temp, CE);
	    strcpy(methodStr, "   (Seaton IP)\n");
	    for (k = 0;  k < Ntemp;  k++)
	      CE[k] *= C3;
	  } else {

	    /* --- Note: Landolt & Boernstein (eq. 34, p. 99) has the
                   wrong power x assigned to neutral and ionic stages,
                   respectively. The proper values should be x = 0.68
                   for neutral stages and x = 0.0 for ionic stages,
                   depending on the form of the appropriate Gaunt factor.
                   Thanks to Jo Bruls for pointing this out -- ------ */

	    strcpy(methodStr, "   (van Regemorter)\n");
	    for (k = 0;  k < Ntemp;  k++)
	      CE[k] = C2_atom / SQ(temp[k]) * f *
		pow((KBOLTZMANN * temp[k])/deltaE, 1.68);
	  }
	  for (k = 0;  k < Ntemp;  k++)
	    fprintf(fp_out, "%9.3E%2s", CE[k],
		    (k == Ntemp-1) ? methodStr : "  ");
	} else {

          /* --- Transition within ionic stage --      -------------- */

	  fprintf(fp_out, " OMEGA  %2d %2d    ", j, i);
	  strcpy(methodStr, "   (van Regemorter)\n");

	  for (k = 0;  k < Ntemp;  k++) {
	    Omega[k] = atom->g[i] * C2_ion / (C1 * temp[k]) * f *
	      ((KBOLTZMANN * temp[k]) / deltaE);
	    fprintf(fp_out, "%9.3E%2s",
		    Omega[k], (k == Ntemp-1) ? methodStr : "  ");
	  }
	}
      }
    }
  }
  /* --- Bound-free transitions
         Modified/JB: Wed Apr 12, 15:43:37, 2000 

         Instead of automatically setting the continuum level (ic)
         to be the ground level of the next higher ionization stage, 
         deduce that from the detailed and fixed B-F transtions 
         present in the model. 
         Since the alpha0 is needed to compute CI anyway, this is 
         a logical way to proceed. --                  -------------- */
    
  for (kr = 0;  kr < atom->Ncont;  kr++) {
    continuum = atom->continuum + kr;
    i  = continuum->i;
    ic = continuum->j;
    alpha0 = continuum->alpha0;
    fprintf(fp_out, " CI     %2d %2d    ", ic, i);
    deltaE = atom->E[ic] - atom->E[i];

    switch (atom->stage[i]) {
    case 0:  gbar_i = 0.1;  break;
    case 1:  gbar_i = 0.2;  break;
    default: gbar_i = 0.3;
    }
    for (k = 0;  k < Ntemp;  k++) {
      CI[k] = C0 / temp[k] * alpha0 * gbar_i *
      ((KBOLTZMANN * temp[k]) / deltaE);
      fprintf(fp_out, "%9.3E%2s", CI[k], (k == Ntemp-1) ? "\n" : "  ");
    }
  }

  for (kf = 0;  kf < atom->Nfixed;  kf++) {
    ft = atom->ft + kf;
    i  = ft->i;
    ic = ft->j;
    if (ft->type == FIXED_CONTINUUM){
      alpha0 = ft->strength;
      fprintf(fp_out, " CI     %2d %2d    ", ic, i);
      deltaE = atom->E[ic] - atom->E[i];

      switch (atom->stage[i]) {
      case 0:  gbar_i = 0.1;  break;
      case 1:  gbar_i = 0.2;  break;
      default: gbar_i = 0.3;
      }
      for (k = 0;  k < Ntemp;  k++) {
	CI[k] = C0 / temp[k] * alpha0 * gbar_i *
	  ((KBOLTZMANN * temp[k]) / deltaE);
	fprintf(fp_out, "%9.3E%2s", CI[k], (k == Ntemp-1) ? "\n" : "  ");
      }
    }
  }
  fprintf(fp_out, "\n END\n");
  free(CE);
}
/* ------- end ---------------------------- crossSection.c ---------- */

/* ------- begin -------------------------- impactParam.c ----------- */

void impactParam(AtomicLine *line, bool_t flag,
		 int Ntemp, double *temp, double *CE)
{
  register int k;

  /* --- Use Gauss-Laguerre quadrature for integration
         over Maxwellian energy distribution of exciting electrons -- */

  int    i = line->i, j = line->j;
  double Ri, Rj, n_eff_min;
  Atom *atom = line->atom;

  /* --- Determine the radius of the orbital --   ------------------- */

  if (flag) {
    n_eff_min = MIN(n_eff[i], n_eff[j]);
    R0 = 0.25 * (5.0*square(n_eff_min) + n_eff_min + 1.0);
  } else {
    Ri = 0.5 * (3.0*square(n_eff[i]) - l[i]*(l[i] + 1.0));
    Rj = 0.5 * (3.0*square(n_eff[j]) - l[j]*(l[j] + 1.0));
    R0 = MIN(Ri, Rj);
  }
  if (VERBOSE >= 1) {
    fprintf(stderr, "R0(%2d, %2d) = %f * a0 (FLAG=%1d)\n",
	    i, j, R0, R_AVG_FLAG);
  }
  /* --- For each temperature integrate over the Maxwellian of the 
         electrons to get the collisional rate coefficient from the
         cross-section in impact parameter approximation -- --------- */

  for (k = 0;  k < Ntemp;  k++) {
    x0    = deltaE / (KBOLTZMANN * temp[k]);
    ERkT  = E_Rydberg / (KBOLTZMANN * temp[k]);
    CE[k] = GaussLaguerre(phiImpact);
  }
}
/* ------- end ---------------------------- impactParam.c ----------- */

/* ------- begin -------------------------- findbeta1.c ------------- */

#define N_MAX_TRY    20
#define DELTA        1.0E-3
#define BETA_START   1.0E-5

double findbeta1(double righthand, double *zeta1)
{

/*  Find solution beta1 of transcendental equation y(beta) = 0.
 *
 *  beta1 -- solution to y(beta) = K0^2(beta) + k1^2(beta) - righthand = 0
 *  zeta1 -- zeta(beta1) = beta1^2 * (K0(beta1)^2 + K1(beta1)^2)
 *                       = beta1^2 * righthand
 *  DELTA -- desired relative accuracy in beta1
 *
 *  Han Uitenbroek, Jan 1992
 */

  int    Ntry;
  double beta, beta1, beta2, y, y_old, delta;

  beta1 = BETA_START;
  y = square(Bessel_K0(beta1)) + square(Bessel_K1(beta1)) - righthand;

  if (y <= 0.0) {
    fprintf(stderr, "Warning: beta1 is smaller than BETA_START\n");
    *zeta1 = righthand * SQ(beta1);
    return beta1;
  }
  /* --- First find turning point by progressively increasing the
         interval [beta1, beta] --                ------------------- */

  Ntry = 0;
  while (Ntry <= N_MAX_TRY) {
    beta2 = 2.0 * beta1;
    y = square(Bessel_K0(beta2)) + square(Bessel_K1(beta2)) - righthand;
    if (y > 0.0) {
      beta1 = beta2;
      Ntry++;
    } else
      break;
  }
  if (Ntry > N_MAX_TRY) {
    fprintf(stderr, "Could not bracket beta1 in %d tries\n", N_MAX_TRY);
    *zeta1 = SQ(beta2) * righthand;
    return beta2;
  }
  /* --- Then narrow solution down by halving the interval until
         relative change in beta1 is less than or equal delta. -- --- */

  while ((delta = 1.0 - beta1/beta2) > DELTA) {
    beta = 0.5*(beta1 + beta2);
    if ((y = square(Bessel_K0(beta)) + square(Bessel_K1(beta))
	 - righthand) > 0.0)
      beta1 = beta;
    else
      beta2 = beta;
  }
  if (VERBOSE >= 2) {
    fprintf(stderr, "Search for beta1 converged with delta = %E\n", delta);
  }
  *zeta1 = square(beta) * righthand;
  return beta;
}
/* ------- end ---------------------------- findbeta1.c ------------- */

/* ------- begin -------------------------- phiImpact.c ------------- */

double phiImpact(double x)
{
  double beta0, beta1, zeta1, sigma_0, sigma_1, righthand;

  x += x0;

  beta0 = sqrt(x/ERkT) * x0/(2*x + x0) * R0;
  sigma_0 = beta0 * Bessel_K0(beta0) * Bessel_K1(beta0);

  righthand = SQ(2*x + x0) / (8 * ERkT * x0 * f);
  beta1 = findbeta1(righthand, &zeta1);
  sigma_1 = 0.5*zeta1 + beta1 * Bessel_K0(beta1) * Bessel_K1(beta1);

  if (VERBOSE >= 1) {
    fprintf(stderr, "sigma_0: %E,  sigma_1: %E\n", sigma_0, sigma_1);
  }
  return 8.0*PIa0square * SQ(ERkT) * (f/x0) * MIN(sigma_0, sigma_1);
}
/* ------- end ---------------------------- phiImpact.c ------------- */

