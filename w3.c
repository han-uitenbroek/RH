/* ------- file: -------------------------- w3.c --------------------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Fri Jan  4 13:30:36 2013 --

       --------------------------                      ----------RH-- */

/* --- Functions for piecewise quadratic integration -- ------------- */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "rh.h"

/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */


/* ------- begin -------------------------- w2.c -------------------- */

void w2(double dtau, double *w)
{
  double expdt;

  if (dtau < 5.0E-4) {
    w[0] = dtau*(1.0 - 0.5*dtau);
    w[1] = SQ(dtau) * (0.5 - dtau/3.0);
  } else if (dtau > 50.0) {
    w[1] = w[0] = 1.0;
  } else {
    expdt = exp(-dtau);
    w[0]  = 1.0 - expdt;
    w[1]  = w[0] - dtau*expdt;
  }
}
/* ------- end ---------------------------- w2.c -------------------- */

/* ------- begin -------------------------- w3.c -------------------- */

void w3(double dtau, double *w)
{
  double expdt, delta;

  if (dtau < 5.0E-4) {
    w[0]   = dtau*(1.0 - 0.5*dtau);
    delta  = SQ(dtau);
    w[1]   = delta*(0.5 - dtau/3.0);
    delta *= dtau;
    w[2]   = delta*(1.0/3.0 - 0.25*dtau);
  } else if (dtau > 50.0) {
    w[1] = w[0] = 1.0;
    w[2] = 2.0;
  } else {
    expdt = exp(-dtau);
    w[0]  = 1.0 - expdt;
    w[1]  = w[0] - dtau*expdt;
    w[2]  = 2.0*w[1] - SQ(dtau) * expdt;
  }
}
/* ------- end ---------------------------- w3.c -------------------- */

/* ------- begin -------------------------- U3.c -------------------- */

void U3(double dtau, double *U)
{
  double expdt, delta;

  if (dtau < 2.0E-4) {
    U[0]   = dtau*(1.0 - 0.5*dtau);
    delta  = SQ(dtau);
    U[1]   = 0.5*delta * (1.0 - dtau/3.0);
    delta *= dtau;
    U[2]   = (delta/3.0) * (1.0 - 0.25*dtau);
  } else if (dtau > 200.0) {
    U[0] = 1.0;
    U[1] = dtau - U[0];
    U[2] = SQ(dtau)  - 2.0*U[1];
  } else {
    U[0]  = 1.0 - exp(-dtau);
    U[1]  = dtau - U[0];
    U[2]  = SQ(dtau) - 2.0*U[1];
  }
}
/* ------- end ---------------------------- U3.c -------------------- */
