/* ------- file: -------------------------- gaussleg.c --------------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Thu Nov  2 10:24:25 2000 --

       --------------------------                      ----------RH-- */

/* --- Compute abscissas and weights for N-point Gauss Legendre
       quadrature formula

  See: Num Rec p. 125 --                               -------------- */

#include <math.h>


#include "rh.h"
#include "constant.h"

/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */


/* ------- begin -------------------------- GaussLeg.c -------------- */

#define ERROR 3.0E-14

void GaussLeg(double x1, double x2, double *x, double *w, int n)
{
  register int j;

  int    i, m;
  double xm, xl, zz, z1, dz, p1, p2, p3, pp;

  m  = (n + 1) / 2;
  xm = 0.5 * (x2 + x1);
  xl = 0.5 * (x2 - x1);

  for (i=0; i<m; i++) {
    zz = cos(PI * (i + 0.75) / (n + 0.5));

    do {
      p1 = 1.0;
      p2 = 0.0;
      for (j=1; j<=n; j++) {
	p3 = p2;
	p2 = p1;
	p1 = (2.0*(j - 0.5)*zz*p2 - (j - 1.0)*p3) / j;
      }
      pp = n * (zz*p1 - p2) / (zz*zz - 1.0);
      z1 = zz;
      zz = z1 - p1/pp;
      dz = fabs(p1/pp);
    } while (dz > ERROR);

    x[i] = (xm - xl*zz);
    w[i] = (2.0 * xl / ((1.0 - zz*zz)*pp*pp));
    x[n-1 - i] = (xm + xl*zz);
    w[n-1 - i] = w[i];
  }
}
/* ------ end ------------------------------ GaussLeg.c ------------- */
