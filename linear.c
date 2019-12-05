/* ------- file: -------------------------- linear.c ----------------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Thu Jun 14 09:53:15 2001 --

       --------------------------                      ----------RH-- */

/* --- Linear interpolation routine. --                -------------- */


#include "rh.h"

/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */


/* ------- begin -------------------------- Linear.c ---------------- */

void Linear(int Ntable, double *xtable, double *ytable,
	    int N, double *x, double *y, bool_t hunt)
{
  register int n;

  bool_t ascend;
  int    j = 0;
  double xmin, xmax, fx;

  ascend = (xtable[1] > xtable[0]) ? TRUE : FALSE;
  xmin = (ascend) ? xtable[0] : xtable[Ntable-1];
  xmax = (ascend) ? xtable[Ntable-1] : xtable[0];

  for (n = 0;  n < N;  n++) {
    if (x[n] <= xmin)
      y[n] = (ascend) ? ytable[0] : ytable[Ntable-1];
    else if (x[n] >= xmax)
      y[n] = (ascend) ? ytable[Ntable-1] : ytable[0];
    else {
      if (hunt) 
	Hunt(Ntable, xtable, x[n], &j);
      else
	Locate(Ntable, xtable, x[n], &j);

      fx = (xtable[j+1] - x[n]) / (xtable[j+1] - xtable[j]);
      y[n] = fx*ytable[j] + (1 - fx)*ytable[j+1];
    }
  }
}
/* ------- end ---------------------------- Linear.c ---------------- */

/* ------- begin -------------------------- BiLinear.c -------------- */

double BiLinear(int Na, double *a_table, double a,
		int Nb, double *b_table, double b,
		double **f, bool_t hunt)
{
  static int i = 0, j = 0;

  double fa, fb;

  /* --- Bi-linear interpolation of function f[][] given on 
         rectangular grid (a_table, b_table) --        -------------- */

  if (hunt) {
    Hunt(Na, a_table, a, &i);
    Hunt(Nb, b_table, b, &j);
  } else {
    Locate(Na, a_table, a, &i);
    Locate(Nb, b_table, b, &j);
  }
  fa = (a_table[i+1] - a) / (a_table[i+1] - a_table[i]);
  fb = (b_table[j+1] - b) / (b_table[j+1] - b_table[j]);

  return                 fa*fb * f[i][j] +
                 fa*(1.0 - fb) * f[i][j+1] +
                 (1.0 - fa)*fb * f[i+1][j] +
         (1.0 - fa)*(1.0 - fb) * f[i+1][j+1];
}
/* ------- end ---------------------------- BiLinear.c -------------- */
