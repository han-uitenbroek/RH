/* ------- file: -------------------------- make_h.c ---------------- */

/* Create Hydrogen atom.

     Cj->i = Ne CE gi/gj sqrt(T)

 * Units: [CE] = s^-1 K^-1/2 m^3
 *
 * Han Uitenbroek
 * Last modified: Thu Feb  3 11:48:44 2000 --
 */
 
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "constant.h"
#include "error.h"
#include "inputs.h"


#define  COMMENT_CHAR  "#"

#define  NTEMP       6
#define  N_MAX_LEVEL 11

#define  H_ABUNDANCE   1.0
#define  H_WEIGHT      1.00790

#define  QCORE_DEFAULT  1.0
#define  QWING_DEFAULT 30.0
#define  NLAMB_DEFAULT   20


/* --- Function prototypes --                          -------------- */

double E1(double x);
double E2(double x);
double f(double x);
double g0(int n);
double g1(int n);
double g2(int n);
double GauntBF(int n);
void   Hboundfree(Atom *H);
void   Johnson_f(Atom *H);
void   Johnson_CE(int Ntemp, double *temp, int Nlevel);
void   Johnson_CI(int Ntemp, double *temp, int Nlevel);
void   writeModelAtom(Atom *atom, FILE *fpOut);


/* --- Global variables --                             -------------- */

double  E_Rydberg;
CommandLine commandline;
char   messageStr[MAX_LINE_SIZE];


/* ------- begin -------------------------- make_H.c ---------------- */

int main(int argc, char *argv[])
{
  register int i;

  char labels[N_MAX_LEVEL][ATOM_LABEL_WIDTH+1] =
    {"H I 1S 2SE          ",
     "H I 2P 2PO          ",
     "H I 3D 2DE          ",
     "H I 4F 2FO          ",
     "H I 5G 2GE          ",
     "H I 6H 2HO          ",
     "H I 7I 2IE          ",
     "H I 8K 2KO          ",
     "H I 9L 2LE          ",
     "H I 10M 2MO         ",
     "H II                " };

  int    Ntemp = NTEMP, Nlevel;
  double temp[NTEMP] = {3.0E+3, 5.0E+3, 7.0E+3, 1.0E+4, 2.0E+4, 3.0E+4};
  Atom H;

  commandline.quiet = FALSE;
  commandline.logfile = stderr;

  if (argc < 2) {
    fprintf(stderr, "Usage: %s Nlevel\n", argv[0]);
    exit(0);
  } else
    Nlevel = atoi(argv[1]);

  if (Nlevel > N_MAX_LEVEL) {
    fprintf(stderr, "Maximum number of levels is %d\n", N_MAX_LEVEL);
    exit(-1);
  }

  strcpy(H.ID, "H ");
  H.abundance = H_ABUNDANCE;
  H.weight    = H_WEIGHT;

  E_Rydberg = E_RYDBERG / (1.0 + M_ELECTRON / (H.weight * AMU));

  H.Nlevel = Nlevel;
  H.Nline  = ((Nlevel - 1) * (Nlevel - 2)) / 2;
  H.Ncont  = Nlevel - 1;
  H.Nfixed = 0;

  H.g = (double *) malloc(H.Nlevel * sizeof(double));
  H.E = (double *) malloc(H.Nlevel * sizeof(double));
  H.stage = (int *) malloc(H.Nlevel * sizeof(int));
  H.label = (char **) matrix_char(H.Nlevel, ATOM_LABEL_WIDTH+1);

  H.line = (AtomicLine *) malloc(H.Nline * sizeof(AtomicLine));
  H.continuum = (AtomicContinuum *) malloc(H.Ncont * sizeof(AtomicContinuum));

  for (i = 0;  i < H.Nlevel-1;  i++) {
    strcpy(H.label[i], labels[i]);
    H.g[i] = 2 * SQ(i+1);
    H.E[i] = E_Rydberg * (1.0 - 1.0/SQ((double) (i+1)));
    H.stage[i] = 0;
  }
  strcpy(H.label[H.Nlevel-1], labels[H.Nlevel-1]);
  H.g[H.Nlevel-1] = 1.0;  H.E[H.Nlevel-1] = E_Rydberg;
  H.stage[H.Nlevel-1] = 1;

  Johnson_f(&H);
  Hboundfree(&H);
  writeModelAtom(&H, stdout);

  Johnson_CE(Ntemp, temp, H.Nlevel);
  Johnson_CI(Ntemp, temp, H.Nlevel);
}
/* ------- end ---------------------------- make_H.c ---------------- */

/* ------- begin -------------------------- Johnson_f.c ------------- */

void Johnson_f(Atom *H)
{
/* Bound-bound oscillator strengths for Hydrogen transitions between
 * states with principal quantum number n and np. 
 *
 * Reference:
 *      -- L.C. Johnson (1972), ApJ 174, 227-236
 */

  register int i, j;

  int    n, np;
  double fnnp, x, C0, C1, lambda0;
  AtomicLine *line = H->line;

  C0 = 2.0*PI * (Q_ELECTRON/EPSILON_0) * (Q_ELECTRON/M_ELECTRON) / CLIGHT;
  C1 = 32.0 / (3.0 * sqrt(3.0) * PI);

  for (i = 0;  i < H->Nlevel-1;  i++) {
    n = i+1;
    for (j = i+1;  j < H->Nlevel-1;  j++) {
      np = j+1;
      x  = 1.0 - SQ((double) n)/SQ((double) np);
      fnnp = C1 * n/CUBE(np*x) * (g0(n) + (g1(n) + g2(n)/x)/x);

      line->i = i;  line->j = j;

      lambda0 = (HPLANCK * CLIGHT) / (H->E[j] - H->E[i]);
      line->Aji = C0 / SQ(lambda0) * (H->g[i] / H->g[j]) * fnnp;
      line->lambda0 = lambda0 / NM_TO_M;
      line->Nlambda = NLAMB_DEFAULT;

      line->symmetric = TRUE;
      line->Voigt = TRUE;
      line->PRD = FALSE;

      line->qcore = QCORE_DEFAULT;
      line->qwing = QWING_DEFAULT;

      line->vdWaals = UNSOLD;
      line->cvdWaals[0] = line->cvdWaals[2] = 1.0;
      line->cvdWaals[1] = line->cvdWaals[3] = 0.0;
      line->Grad = line->cStark = 0.0;

      line++;
    }
  }
}
/* ------- end ---------------------------- Johnson_f.c ------------- */

/* ------- begin -------------------------- Hboundfree.c ------------ */

void Hboundfree(Atom *H)
{
  register int i;

  double lambda0, C0;
  AtomicContinuum *continuum = H->continuum;

  C0 = 32.0/(3.0*sqrt(3.0)) * SQ(Q_ELECTRON/sqrt(4.0*PI*EPSILON_0)) /
    (M_ELECTRON * CLIGHT) * HPLANCK / (2.0*E_Rydberg);

  for (i = 0;  i < H->Nlevel-1;  i++) {
    continuum->hydrogenic = TRUE;

    continuum->i = i;  continuum->j = H->Nlevel-1;
    lambda0 = (HPLANCK * CLIGHT) * SQ((double) (i+1)) / E_Rydberg;
    continuum->lambda = (double *) malloc(continuum->Nlambda * sizeof(double));
    continuum->lambda[0] = 0.25*lambda0 / NM_TO_M;

    continuum->alpha0 = C0 * (i+1) * GauntBF(i+1);
    continuum->Nlambda = NLAMB_DEFAULT;
 
    continuum++;
  }
}
/* ------- end ---------------------------- Hboundfree.c ------------ */

/* ------- begin -------------------------- Johnson_CI.c ------------ */

/* Collisional ionization rate coefficients for state with principal
 * quantum number n.
 * 
 * Define CI through:  Cn->cont = Ne * CI * exp(-dE/kT) * sqrt(T)
 *
 * Reference:
 *      -- L.C. Johnson (1972), ApJ 174, 227-236
 */

void Johnson_CI(int Ntemp, double *temp, int Nlevel)
{
  register int k, i;

  int    n;
  double C0, C1, PIa0sq, rn, y, z, kT_m, Gn, Bn, bn, CI;

  C0 = sqrt((8 * KBOLTZMANN) / (PI * M_ELECTRON));
  C1 = 32.0 / (3.0 * sqrt(3.0) * PI);
  PIa0sq = PI * SQ(RBOHR);

  /* --- Print used temperatures --                    -------------- */

  fprintf(stdout, "\n%1s Collisional rate coefficients from Johnson_CI.c\n\n",
	  COMMENT_CHAR);
  fprintf(stdout, " TEMP   %2d       ", Ntemp);
  for (k = 0;  k < Ntemp;  k++)
    fprintf(stdout, "%9.1f%2s", temp[k], (k == Ntemp-1) ? "\n\n" : "  ");

  for (i = 0;  i < Nlevel - 1;  i++) {
    n = i + 1;

    if (n == 1) {
      rn =  0.45;
      bn = -0.603;
    } else {
      rn = 1.94*pow(n, -1.57);
      bn = (4.0 + (-18.63 + (36.24 - 28.09/n)/n)/n) / n;
    }
    Gn = C1*n * (g0(n)/3.0 + g1(n)/4.0 + g2(n)/5.0); 
    Bn = 2.0*SQ(n)/3.0 * (5.0 + bn);

    fprintf(stdout, " CI     %2d %2d    ", i, Nlevel - 1);
    for (k = 0;  k < Ntemp;  k++) {
      kT_m = C0 * sqrt(temp[k]);
      y = E_Rydberg / (n*n * KBOLTZMANN *temp[k]);
      z = rn + y;

      CI  = 2.0*PIa0sq * SQ(n) * kT_m * SQ(y)*(Gn*(E1(y)/y - E1(z)/z) +
		      (Bn - Gn*log(2.0*SQ(n))) * (f(y) - f(z)));
      CI *= exp(y) / sqrt(temp[k]);

      fprintf(stdout, "%9.3e%2s", CI,
	      (k == Ntemp-1) ? "   (Johnson)\n" : "  ");
    }
  }
}
/* ------- end ---------------------------- Johnson_CI.c ------------ */

/* ------- begin -------------------------- Johnson_CE.c ------------ */

/* Collisional excitation rate coefficients for state with principal
 * quantum number n to np.
 * 
 * Define CI through:  Cn->np = Ne * CE * exp(-dE/kT) * sqrt(T)
 *
 * Reference:
 *      -- L.C. Johnson (1972), ApJ 174, 227-236
 */

void Johnson_CE(int Ntemp, double *temp, int Nlevel)
{
  register int k, i, j;

  int    n, np;
  double C0, C1, PIa0sq, x, y, z, rnnp, kT_m, Annp, fnnp, Bnnp, bn, CE;

  C0 = sqrt((8 * KBOLTZMANN) / (PI * M_ELECTRON));
  C1 = 32.0 / (3.0 * sqrt(3.0) * PI);
  PIa0sq = PI * SQ(RBOHR);

  /* --- Print used temperatures --                    -------------- */

  fprintf(stdout, "\n%1s Collisional rate coefficients from Johnson_CE.c\n\n",
	  COMMENT_CHAR);
  fprintf(stdout, " TEMP   %2d       ", Ntemp);
  for (k = 0;  k < Ntemp;  k++)
    fprintf(stdout, "%9.1f%2s", temp[k], (k == Ntemp-1) ? "\n\n" : "  ");

  for (i = 0;  i < Nlevel - 1;  i++) {
    n = i + 1;
    if (n == 1) {
      rnnp =  0.45;
      bn   = -0.603;
    } else {
      rnnp = 1.94*pow(n, -1.57);
      bn   = (4.0 + (-18.63 + (36.24 - 28.09/n)/n)/n) / n;
    }
    for (j = i + 1;  j < Nlevel - 1;  j++) {
      np = j + 1;
      x  = 1.0 - SQ((double) n)/SQ((double) np);

      rnnp *= x;
      fnnp = C1 * n/CUBE(np*x) * (g0(n) + (g1(n) + g2(n)/x)/x);
      Annp = 2.0*SQ(n)/x * fnnp;
      Bnnp = 4.0*SQ(SQ(n))/(CUBE(np) * SQ(x)) *
	(1.0 + 4.0/(3.0*x) + bn/SQ(x));

      fprintf(stdout, " CE     %2d %2d    ", j, i);
      for (k = 0;  k < Ntemp;  k++) {
	kT_m = C0 * sqrt(temp[k]);
	y = x * E_Rydberg / (SQ(n) * KBOLTZMANN * temp[k]);
	z = rnnp + y;

	CE = 2*PIa0sq * SQ(n)/x * kT_m * SQ(y) *
	  (Annp * (E1(y)*(1.0/y + 0.5) - E1(z)*(1.0/z + 0.5)) +
	   (Bnnp - Annp*log(2.0*SQ(n)/x)) * (E2(y)/y - E2(z)/z));
	CE *= exp(y) / sqrt(temp[k]);

	fprintf(stdout, "%9.3e%2s", CE,
		(k == Ntemp-1) ? "   (Johnson)\n" : "  ");
      }
    }
  }
}
/* ------- end ---------------------------- Johnson_CI.c ------------ */

/* ------- begin -------------------------- f.c --------------------- */

double f(double x)
{
  return exp(-x)/x - 2*E1(x) + E2(x);
}
/* ------- end ---------------------------- f.c --------------------- */

/* ------- begin -------------------------- g0.c -------------------- */

double g0(int n)
{
  switch (n) {
  case 1:  return 1.133;
  case 2:  return 1.0785;
  default: return 0.9935 + (0.2328 - 0.1296/n)/n;
  }
}
/* ------- end ---------------------------- g0.c -------------------- */

/* ------- begin -------------------------- g1.c -------------------- */

double g1(int n)
{
  switch (n) {
  case 1:  return -0.4059;
  case 2:  return -0.2319;
  default: return -(0.6282 - (0.5598 - 0.5299/n)/n) / n;
  }
}
/* ------- end ---------------------------- g1.c -------------------- */

/* ------- begin -------------------------- g2.c -------------------- */

double g2(int n)
{
  switch (n) {
  case 1:  return  0.07014;
  case 2:  return  0.02947;
  default: return (0.3887 - (1.181 - 1.4700/n)/n) / (n*n);
  }
}
/* ------- end ---------------------------- g2.c -------------------- */

/* ------- begin -------------------------- GauntBF.c --------------- */

double GauntBF(int n)
{
  /* --- M. J. Seaton (1960), Rep. Prog. Phys. 23, 313 -- ----------- */

  double x, x3, nsqx;

  x    = 1.0 / SQ((double) n);
  x3   = pow(x, 0.33333333);
  nsqx = 1.0 / (SQ(n) * x);

  return 1.0 + 0.1728*x3 * (1.0 - 2.0*nsqx) -
               0.0496*SQ(x3) * (1.0 - (1.0 - nsqx)*0.66666667*nsqx);
}
/* ------- end ---------------------------- GauntBF.c --------------- */

