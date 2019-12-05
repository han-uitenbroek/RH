/* ------- file: -------------------------- hcollisions.c ----------- */

/* Collisional excitation and ionization routine for hydrogen.
 * Prints out CE, where CE is defined through:

     Cj->i = Ne CE gi/gj sqrt(T)

 * Units: [CE] = s^-1 K^-1/2 m^3
 *
 * Han Uitenbroek
 * Last modified: Thu Feb  3 11:41:55 2000 --
 */
 
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "rh.h"
#include "constant.h"
#include "error.h"
#include "inputs.h"

#define  NTEMP  6


/* --- Function prototypes --                          -------------- */

double crossSection(double x);
int    get_ilevel(char *state);
double GaussLaguerre(double (*function) (double x));
void   Giovanardi(int Ntemp, double *temp);
void   oneS_twoS(int Ntemp, double *temp);
void   oneS_twoP(int Ntemp, double *temp);
void   Whelan(int Ntemp, double *temp);


/* --- Global variables --                             -------------- */

CommandLine commandline;
char   messageStr[MAX_LINE_SIZE];


/* ------- begin -------------------------- Hcollisions.c ----------- */

int main(int argc, char *argv[])
{
  register int k, n;

  int    Ntemp = NTEMP;
  double temp[NTEMP] = {3.0E+3, 5.0E+3, 7.0E+3, 1.0E+4, 2.0E+4, 3.0E+4};

  commandline.quiet = FALSE;
  commandline.logfile = stderr;

  /* --- Print used temperatures --                    -------------- */

  fprintf(stdout, "\n* Collisional rate coefficients from hcollisions.c\n\n");
  fprintf(stdout, " TEMP   %2d       ", Ntemp);
  for (k = 0;  k < Ntemp;  k++)
    fprintf(stdout, "%9.1f%2s", temp[k], (k == Ntemp-1) ? "\n\n" : "  ");

  oneS_twoS(Ntemp, temp);
  oneS_twoP(Ntemp, temp);
  Giovanardi(Ntemp, temp);
  Whelan(Ntemp, temp);
}
/* ------- end ---------------------------- Hcollisions.c ----------- */

/* ------- begin -------------------------- oneS_twoS.c ------------- */

/* Effective collision strength for Hydrogen 1s - 2s.
 *
 * Reference:
 *   -- T.T. Scholz, H.R.J. Walters, P.G. Burke, and M.P. Scott (1990),
 *      MNRAS 242, 692-697
 */

#define NB   8

void oneS_twoS(int Ntemp, double *temp)
{
  register int k;

  static double b[NB] = { 4.5168E-02,  2.8056E+01,  7.2945E+00,
                          2.4805E-01,  1.0044E-01, -1.1143E-02,
                         -1.3432E-03,  3.7570E-04 };
  double x, gi = 2.0, C1, CE;

  C1 = ((E_RYDBERG/sqrt(M_ELECTRON)) * PI*SQ(RBOHR)) *
    sqrt(8.0/(PI*KBOLTZMANN));
  fprintf(stdout, " CE     %2d %2d    ", 1, 0);

  for (k = 0;  k < Ntemp;  k++) {
    x = (KBOLTZMANN * temp[k]) / E_RYDBERG;
    CE  = b[0]*log(b[1] * x)*exp(-b[2]*x) +
          b[3] + x*(b[4] + x*(b[5] + x*(b[6] + x*b[7])));
    CE *= C1 / (gi * temp[k]);
    fprintf(stdout, "%9.3e%2s", CE,
	    (k == Ntemp-1) ? "   (Scholz et al.)\n" : "  ");
  }
}
/* ------- end ---------------------------- oneS_twoS.c ------------- */

/* ------- begin -------------------------- oneS_twoP.c ------------- */

/* Effective collision strength for Hydrogen 1s - 2p.
 *
 * Reference:
 *   -- T.T. Scholz, H.R.J. Walters, P.G. Burke, and M.P. Scott (1990),
 *      MNRAS 242, 692-697
 */

#define NC   6

void oneS_twoP(int Ntemp, double *temp)
{
  register int k;

  static double c[NC] = {  3.6177E-01,  1.3891E+00,  5.0866E-01,
			  -3.8011E-01,  1.0158E-01, -1.0072E-02 };
  double x, gi = 2.0, C1, CE;

  C1 = ((E_RYDBERG/sqrt(M_ELECTRON)) * PI*SQ(RBOHR)) *
    sqrt(8.0/(PI*KBOLTZMANN));
  fprintf(stdout, " CE     %2d %2d    ", 2, 0);

  for (k = 0;  k < Ntemp;  k++) {
    x = (KBOLTZMANN * temp[k]) / E_RYDBERG;
    CE  = c[0] + x*(c[1] + x*(c[2] + x*(c[3] + x*(c[4] + x*c[5]))));
    CE *= C1 / (gi * temp[k]);
    fprintf(stdout, "%9.3e%2s", CE,
	    (k == Ntemp-1) ? "   (Scholz et al.)\n" : "  ");
  }
}
/* ------- end ---------------------------- oneS_twoP.c ------------- */

/* ------- begin -------------------------- Giovanardi.c ------------ */

/* Electron impact excitation rates for Hydrogen.
 *
 * References:
 *      -- C. Giovanardi, A. Natta, and F.Palla (1987),
 *         Astron. Astrophys. Suppl. 70, 269-280
 *      -- C. Giovanardi, F. Palla (1989), 
 *         Astron. Astrophys. Suppl. 77, 157-160
 *
 * See also:
 *         E.S Chang, E.H. Avrett, and R. Loeser (1991), A&A 247, 580-583
 */

#define  GIOVANARDI_DATA  "/home/uitenbr/src/rh/Atoms/Giovanardi_etal.dat"

#define  NDATA     21
#define  NLEVEL    11
#define  NPOLYNOME  4

#define  TMIN   2.0E+03
#define  T1     6.0E+03
#define  T2     5.5E+04
#define  T3     7.2E+04
#define  TMAX   5.0E+05

void Giovanardi(int Ntemp, double *temp)
{
  register int k;

  char   line1[MAX_LINE_SIZE], line2[MAX_LINE_SIZE];
  bool_t exit_on_EOF;
  int    i, j, ij, Nlevel = NLEVEL;
  double a[NPOLYNOME], b[NPOLYNOME], corr, T, fraction, 
         g[NLEVEL] = {2.0, 2.0, 6.0, 2.0, 6.0, 10.0,
                      2.0, 6.0, 10.0, 14.0, 50.0}, Clo, Chi, C1, **CE;
  FILE *fp_data;

  if ((fp_data = fopen(GIOVANARDI_DATA, "r")) == NULL) {
    sprintf(messageStr, "Unable to open input file %s", GIOVANARDI_DATA);
    Error(ERROR_LEVEL_2, "Giovanardi", messageStr);
  }

  CE = (double **) matrix_double(Nlevel*Nlevel, Ntemp);

  while (getLine(fp_data, "*", line1, exit_on_EOF=FALSE) != EOF) {
    i = get_ilevel(strtok(line1, "-"));
    j = get_ilevel(strtok(NULL, "\n"));

    getLine(fp_data, "*", line1, exit_on_EOF=TRUE);
    getLine(fp_data, "*", line2, exit_on_EOF=TRUE);

    if ((i != -1)  &&  (j != -1)) {
      ij = i*NLEVEL + j;
      sscanf(line1, "%lf %lf %lf %lf %lf", a, a+1, a+2, a+3, &corr);
      sscanf(line2, "%lf %lf %lf %lf", b, b+1, b+2, b+3);

      for (k = 0;  k < Ntemp;  k++) {
        T = temp[k];
	Clo = a[0] + T*(a[1] + T*(a[2] + T*a[3]));
	Chi = b[0] + T*(b[1] + T*(b[2] + T*b[3]));

	if ((T >= T3)  &&  (T <= TMAX)) 
	  CE[ij][k] += Chi;
	else if ((T > T2)  &&  (T < T3)) {
	  fraction = (T3 - T) / (T3 - T2);
	  CE[ij][k] += fraction*Clo + (1.0 - fraction)*Chi;
	} else if ((T >= T1)  &&  (T <= T2))
	  CE[ij][k] += Clo;
	else if ((T >= TMIN)  &&  (T < T1))

	  /* --- Correction for low temperature --     -------------- */

	  CE[ij][k] += Clo * (1.0 + (corr/100.0)*pow(6000.0 - T, 2.2));
      }
    }
  }
  /* --- Now print the results --                      -------------- */

  C1 = ((E_RYDBERG/sqrt(M_ELECTRON)) * PI*SQ(RBOHR)) *
    sqrt(8.0/(PI*KBOLTZMANN));

  for (i = 0;  i < NLEVEL;  i++) {
    for (j = i;  j < NLEVEL;  j++) {
      ij = i*NLEVEL + j;
      if (CE[ij][0] > 0.0) {
	fprintf(stdout, " CE     %2d %2d    ", j, i);
        for (k = 0;  k < Ntemp;  k++) {
	  CE[ij][k] *= C1 / (g[i] * temp[k]);
	  fprintf(stdout, "%9.3e%2s", CE[ij][k],
	    (k == Ntemp-1) ? "   (Giovanardi et al.)\n" : "  ");
	}
      }
    }
  }
  freeMatrix((void **) CE);
  fclose(fp_data);
}
/* ------- end ---------------------------- Giovanardi.c ------------ */

/* ------- begin -------------------------- get_ilevel.c ------------ */

int get_ilevel(char *thestate)
{
  register int i;

  static char state[NDATA][3] = {"1s", "2s", "2p", "3s", "3p", "3d",
                                "4s", "4p", "4d", "4f", "5", "6", "7",
                                "8", "9", "10", "11", "12", "13", "14", "15"};
  static int ilevel[NDATA] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
                              -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};

  for (i = 0;  i < NDATA;  i++)
    if (strstr(state[i], thestate)) break;

  return ilevel[i];
}
/* ------- end ---------------------------- get_ilevel.c ------------ */

/* ------- begin -------------------------- Whelan.c ---------------- */

#define NW              6
#define SIGMA_LOW    10.0
#define SIGMA_HIGH   20.0

#define  WHELAN_DATA  "/home/uitenbr/src/rh/Atoms/Whelan.dat"

double x0;

void Whelan(int Ntemp, double *temp)
{
  register int k, l, n;

  char   line[MAX_LINE_SIZE], labeli[3], labelj[3];
  bool_t exit_on_EOF;
  int    i, j, nread, nQ, nQp;
  double  k1[NW] = {1.0, 1.2, 1.606238, 2.0, 2.711088, 3.0}, 
         x[NW], tension, cs[NW], C0, CE;
  FILE  *fp_data;

  C0 = sqrt(8*KBOLTZMANN / (PI*M_ELECTRON)) * PI * SQ(RBOHR);

  if ((fp_data = fopen(WHELAN_DATA, "r")) == NULL) {
    sprintf(messageStr, "Unable to open input file %s", WHELAN_DATA);
    Error(ERROR_LEVEL_2, "Whelan", messageStr);
  }

  getLine(fp_data, "*", line, exit_on_EOF=TRUE);
  sscanf(line, "%lf %lf %lf %lf %lf %lf", k1, k1+1, k1+2, k1+3, k1+4, k1+5);

  while (getLine(fp_data, "*", line, exit_on_EOF=FALSE) != EOF) {
    nread = sscanf(line, "%2s-%2s %lf %lf %lf %lf %lf %lf",
		   labeli, labelj, cs, cs+1, cs+2, cs+3, cs+4, cs+5);

    j = get_ilevel(labelj);
    i = get_ilevel(labeli);

    if ((labeli[0] == '4') || (labeli[0] == '3' && labeli[1] == labelj[1]))
      tension = SIGMA_HIGH;
    else
      tension = SIGMA_LOW;
 
    fprintf(stdout, " CE     %2d %2d    ", j, i);
    nQ  = (int) (labeli[0] - '0');
    nQp = (int) (labelj[0] - '0');

    for (k = 0;  k < Ntemp;  k++) {
      x0 = (1.0/(nQ*nQ) - 1.0/(nQp*nQp)) * E_RYDBERG / (KBOLTZMANN * temp[k]);
      for (n = 0;  n < NW;  n++)
	x[n] = (SQ(k1[n]) * E_RYDBERG) / (KBOLTZMANN * temp[k]);

      exp_splineCoef(NW, x, cs, tension);
      CE = C0 * GaussLaguerre(crossSection);
      fprintf(stdout, "%9.3e%2s", CE, (k == Ntemp-1) ? "   (Whelan)\n" : "  ");
    }
  }
}
/* ------- end ---------------------------- Whelan.c ---------------- */

/* ------- begin -------------------------- crossSection.c ---------- */

double crossSection(double x)
{
  double  xT, cs;

  xT = x + x0;
  exp_splineEval(1, &xT, &cs, FALSE);
  return  xT * cs;
}
/* ------- end ---------------------------- crossSection.c ---------- */
