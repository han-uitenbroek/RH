/* ------- file: -------------------------- expspline.c -------------
       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Tue Feb 16 14:37:25 1999 --

       --------------------------                      ----------RH-- */

/* --- Exponential spline interpolation routine.

  Ref: J. Stoer and R. Bulirsch, in Introduction to Numerical
       Analysis, Springer-Verlag, p. 97-102, and excercise 2.32 on p. 115)
       --                                              -------------- */
 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "rh.h"

/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

static bool_t  ascend;
static int     Ntable;
static double *xtable, xmin, xmax, sigma, *M = NULL,
              *ytable, *sinhh = NULL;


/* ------- begin -------------------------- exp_splineCoef.c -------- */

void exp_splineCoef(int N, double *x, double *y, double tension)
{
  register int j;
  static double *u = NULL;

  double *q, p, h, sh, Aj, Cj, Dj, Dj1, Bj, Bj1;

  ascend = (x[1] > x[0]) ? TRUE : FALSE;
  xmin = (ascend) ? x[0] : x[N-1];
  xmax = (ascend) ? x[N-1] : x[0];

  /* --- Use uniform tension scaled to size of table grid -- -------- */

  sigma = tension * (N - 1)/fabs(x[N-1] - x[0]);

  q = M = (double *) realloc(M, N * sizeof(double));
  u = (double *) realloc(u, N * sizeof(double));

  /* --- Store sinh(sigma*h) for use in evaluation part -- ---------- */

  sinhh = (double *) realloc(sinhh, (N - 1)*sizeof(double));

  h  = (x[1] - x[0]);
  sh = sigma*h;
  sinhh[0] = sinh(sh);
  Aj = (1 - sh/sinhh[0]) / h;
  Bj = (sh*cosh(sh)/sinhh[0] - 1) / h; 
  Dj = (y[1] - y[0]) / h;

  /* --- Gaussian elimination --                       -------------- */

  q[0] = u[0] = 0.0;
  for (j = 1;  j < N-1;  j++) {
    h  = x[j+1] - x[j];
    sh = sigma * h;
    sinhh[j] = sinh(sh);
    Bj1 = (sh*cosh(sh)/sinhh[j] - 1) / h;
    Cj  = (1 - sh/sinhh[j]) / h;
    Dj1 = (y[j+1] - y[j]) / (x[j+1] - x[j]);

    p = Aj*q[j-1] + Bj + Bj1;
    q[j] = -Cj / p;
    u[j] = (Dj1 - Dj - Aj*u[j-1]) / p;

    Aj = Cj;  Bj = Bj1;  Dj = Dj1;
  }
  /* --- Backsubstitution --                          --------------- */

  M[N - 1] = 0.0;
  for (j = N-2;  j >= 0;  j--) {
    M[j] = q[j]*M[j+1] + u[j];
  }
  /* --- Store table for use in evaluation part --    --------------- */

  Ntable = N;
  xtable = x;  ytable = y;
}
/* ------- end ---------------------------- exp_splineCoef.c -------- */

/* ------- begin -------------------------- exp_splineEval.c -------- */

void exp_splineEval(int N, double *x, double *y, bool_t hunt)
{
  register int n;

  int    j = 0;
  double hj, fx;

  for (n = 0;  n < N;  n++) {

    /* --- Locate the position of x[n] in table --  ----------------- */

    if (x[n] <= xmin)
      y[n] = (ascend) ? ytable[0] : ytable[Ntable-1];
    else if (x[n] >= xmax)
      y[n] = (ascend) ? ytable[Ntable-1] : ytable[0];
    else {
      if (hunt) 
	Hunt(Ntable, xtable, x[n], &j);
      else
	Locate(Ntable, xtable, x[n], &j);

      /* --- Evaluate exponential spline --             ------------- */

      hj = xtable[j+1] - xtable[j];
      fx = (x[n] - xtable[j]) / hj;

      y[n] = (1 - fx)*ytable[j] + fx*ytable[j+1] +
         (sinh(sigma*(xtable[j+1] - x[n]))/sinhh[j] - 1 + fx) * M[j] +
         (sinh(sigma*(x[n] - xtable[j]))/sinhh[j] - fx) * M[j+1];
    }
  }
}
/* ------- end ---------------------------- exp_splineEval.c -------- */
