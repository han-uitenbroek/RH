/* ------- file: -------------------------- anglequad.c -------------

       Version:       rh2.0, 3-D Cartesian, short characteristics
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Tue Apr 21 14:38:18 2009 --

       --------------------------                      ----------RH-- */

/* --- Determine angular quadrature and fill ray stucture AngleQuadr.

  See: B. G. Carlson in "Methods of Computational Physics", Vol 1,
       (Edited by B. Alder, S. Fernbach, and M. Rottenberg).
       Academic Press, New York (1963).

       By setting ANGLESET = SET_GL_NiXNa in keyword.input a quadrature
       can be chosen that has, per octant, Ni Gauss-Legendre quadrature
       points in inclination, and Na equidistant angle points in
       azimuth (providing Ni x Na rays per octant). 

       Weights are normalized to yield 1.0 when integrated over four
       octants, which represent all space.
       --                                              -------------- */

 
#include <stdlib.h>
#include <math.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "geometry.h"
#include "error.h"
#include "constant.h"

/* --- Define number of rays per octant (NRO) for each of the Carlson
       sets --                                         -------------- */

#define NRO_A2   1
#define NRO_A4   3
#define NRO_A6   6
#define NRO_A8   10
#define NRO_B4   3
#define NRO_B6   6
#define NRO_B8   10

#define HALFPI  PI/2.0


/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern char messageStr[];


/* ------- begin -------------------------- getAngleQuadr.c --------- */

void getAngleQuadr(Geometry *geometry)
{
  const char routineName[] = "getAngleQuadr";
  register int mu, n, m;

  double wnorm, *wmu, *xmu, alpha, singamma;

  /* --- Set A2 --                                     -------------- */

  double mux_A2[NRO_A2] = {0.57735026};
  double muy_A2[NRO_A2] = {0.57735026};
  double wmu_A2[NRO_A2] = {1.00000000};

  /* --- Set A4 --                                     -------------- */

  double mux_A4[NRO_A4] = {0.88191710, 0.33333333, 0.33333333};
  double muy_A4[NRO_A4] = {0.33333333, 0.88191710, 0.33333333};
  double wmu_A4[NRO_A4] = {0.33333333, 0.33333333, 0.33333333};

  /* --- Set A6 --                                     -------------- */

  double mux_A6[NRO_A6] = {0.93094934, 0.68313005, 0.25819889,
			   0.68313005, 0.25819889, 0.25819889};
  double muy_A6[NRO_A6] = {0.25819889, 0.68313005, 0.93094934,
			   0.25819889, 0.68313005, 0.25819889};
  double wmu_A6[NRO_A6] = {0.18333333, 0.15000000, 0.18333333,
			   0.15000000, 0.15000000, 0.18333333};

  /* --- Set A8 --                                     -------------- */

  double mux_A8[NRO_A8] = {0.95118973, 0.78679579, 0.57735027,
			   0.21821789, 0.78679579, 0.57735027,
			   0.21821789, 0.57735027, 0.21821789, 0.21821789};
  double muy_A8[NRO_A8] = {0.21821789, 0.57735027, 0.78679579,
			   0.95118973, 0.21821789, 0.57735027,
			   0.78679579, 0.21821789, 0.57735027, 0.21821789};
  double wmu_A8[NRO_A8] = {0.12698138, 0.09138353, 0.09138353,
			   0.12698138, 0.09138353, 0.07075469,
			   0.09138353, 0.09138353, 0.09138353, 0.12698138};

  /* --- Set B4 --                                     -------------- */

  double mux_B4[NRO_B4] = {0.70412415, 0.09175171, 0.09175171};
  double muy_B4[NRO_B4] = {0.09175171, 0.70412415, 0.09175171};
  double wmu_B4[NRO_B4] = {0.33333333, 0.33333333, 0.33333333};

  /* --- Set B6
         Weights from Jo Bruls, Apr 26 1999 --         -------------- */

  double mux_B6[NRO_B6] = {0.80847426, 0.57735027, 0.11417547,
			   0.57735027, 0.11417547, 0.11417547};
  double muy_B6[NRO_B6] = {0.11417547, 0.57735027, 0.80847426,
			   0.11417547, 0.57735027, 0.11417547};
  double wmu_B6[NRO_B6] = {0.22222222, 0.11111111, 0.22222222,
			   0.11111111, 0.11111111, 0.22222222};

  /* --- Set B8 --                                     -------------- */

  double mux_B8[NRO_B8] = {0.85708018, 0.70273364, 0.50307327,
			   0.11104440, 0.70273364, 0.50307327,
			   0.11104440, 0.50307327, 0.11104440, 0.11104440};
  double muy_B8[NRO_B8] = {0.11104440, 0.50307327, 0.70273364,
			   0.85708018, 0.11104440, 0.50307327,
			   0.70273364, 0.11104440, 0.50307327, 0.11104440};
  double wmu_B8[NRO_B8] = {0.17152232, 0.06900520, 0.06900520,
			   0.17152232, 0.06900520, 0.07140185,
			   0.06900520, 0.06900520, 0.06900520,
			   0.17152232};

  switch (atmos.angleSet.set) {
  case SET_VERTICAL: geometry->Nrays = 1;   break;
  case SET_GL: geometry->Nrays =
    4*atmos.angleSet.Ninclination * atmos.angleSet.Nazimuth;
    break;
  case SET_A2: geometry->Nrays = 4*NRO_A2;  break;
  case SET_A4: geometry->Nrays = 4*NRO_A4;  break;
  case SET_A6: geometry->Nrays = 4*NRO_A6;  break;
  case SET_A8: geometry->Nrays = 4*NRO_A8;  break;
  case SET_B4: geometry->Nrays = 4*NRO_B4;  break;
  case SET_B6: geometry->Nrays = 4*NRO_B6;  break;
  case SET_B8: geometry->Nrays = 4*NRO_B8;  break;
  case NO_SET:
    sprintf(messageStr, "No angleset was specified: %d\n", NO_SET);
    Error(ERROR_LEVEL_2, routineName, messageStr);
    break;
  default:
    sprintf(messageStr, "Did not recognize angleset: %d\n",
	    atmos.angleSet.set);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }

  geometry->mux = (double *) malloc(geometry->Nrays * sizeof(double));
  geometry->muy = (double *) malloc(geometry->Nrays * sizeof(double));
  geometry->muz = (double *) malloc(geometry->Nrays * sizeof(double));
  geometry->wmu = (double *) malloc(geometry->Nrays * sizeof(double));

  switch (atmos.angleSet.set) {
  case SET_VERTICAL:
    geometry->mux[0] = 0.0;
    geometry->muy[0] = 0.0;
    geometry->muz[0] = 1.0;
    geometry->wmu[0] = 1.0;

    return;
  case SET_GL:

    /* --- Quadrature in inclination --                -------------- */

    wmu = (double *) malloc(atmos.angleSet.Ninclination *
			    sizeof(double));
    xmu = (double *) malloc(atmos.angleSet.Ninclination *
			    sizeof(double));
    GaussLeg(0.0, 1.0, xmu, wmu, atmos.angleSet.Ninclination);

    /* --- Grid step of quadrature in azimuth --       -------------- */

    alpha = HALFPI / (2.0 * atmos.angleSet.Nazimuth);

    /* --- Fill first octant --                        -------------- */

    for (n = 0, mu = 0;  n < atmos.angleSet.Ninclination;  n++) {
      singamma = sqrt(1.0 - SQ(xmu[n]));
      for (m = 0;  m < atmos.angleSet.Nazimuth;  m++, mu++) {
	geometry->mux[mu] = cos((2.0*m + 1.0) * alpha) * singamma;
	geometry->muy[mu] = sin((2.0*m + 1.0) * alpha) * singamma;
	geometry->wmu[mu] = wmu[n] / atmos.angleSet.Nazimuth;
      }
    }

    free(wmu);
    free(xmu);
    break;
  case SET_A2:
    for (mu = 0;  mu < NRO_A2;  mu++) {
      geometry->mux[mu] = mux_A2[mu];
      geometry->muy[mu] = muy_A2[mu];
      geometry->wmu[mu] = wmu_A2[mu];
    }
    break;
  case SET_A4:
    for (mu = 0;  mu < NRO_A4;  mu++) {
      geometry->mux[mu] = mux_A4[mu];
      geometry->muy[mu] = muy_A4[mu];
      geometry->wmu[mu] = wmu_A4[mu];
    }
    break;
  case SET_A6:
    for (mu = 0;  mu < NRO_A6;  mu++) {
      geometry->mux[mu] = mux_A6[mu];
      geometry->muy[mu] = muy_A6[mu];
      geometry->wmu[mu] = wmu_A6[mu];
    }
    break;
  case SET_A8:
    for (mu = 0;  mu < NRO_A8;  mu++) {
      geometry->mux[mu] = mux_A8[mu];
      geometry->muy[mu] = muy_A8[mu];
      geometry->wmu[mu] = wmu_A8[mu];
    }
    break;
  case SET_B4:
    for (mu = 0;  mu < NRO_B4;  mu++) {
      geometry->mux[mu] = mux_B4[mu];
      geometry->muy[mu] = muy_B4[mu];
      geometry->wmu[mu] = wmu_B4[mu];
    }
    break;
  case SET_B6:
    for (mu = 0;  mu < NRO_B6;  mu++) {
      geometry->mux[mu] = mux_B6[mu];
      geometry->muy[mu] = muy_B6[mu];
      geometry->wmu[mu] = wmu_B6[mu];
    }
    break;
  case SET_B8:
    for (mu = 0;  mu < NRO_B8;  mu++) {
      geometry->mux[mu] = mux_B8[mu];
      geometry->muy[mu] = muy_B8[mu];
      geometry->wmu[mu] = wmu_B8[mu];
    }
    break;
  default:;
  }
  /* --- Normalize, remember we integrate over four octants -- ------ */

  wnorm = 0.0;
  for (mu = 0;  mu < geometry->Nrays/4;  mu++)
    wnorm += geometry->wmu[mu];
  wnorm = 0.25 / wnorm;

  /* --- Fill the second octant with rays reflected in the y-z plane  */

  for (mu = 0;  mu < geometry->Nrays/4;  mu++) {
    geometry->mux[geometry->Nrays/4 + mu]  = -geometry->mux[mu];
    geometry->muy[geometry->Nrays/4 + mu]  =  geometry->muy[mu];
    geometry->wmu[mu]                     *=  wnorm;
    geometry->wmu[geometry->Nrays/4 + mu]  =  geometry->wmu[mu];
  }
  /* --- Fill the third and fourth octant with rays reflected in
         the x-z plane --                              -------------- */

  for (mu = 0;  mu < geometry->Nrays/2;  mu++) {
    geometry->mux[geometry->Nrays/2 + mu]  = -geometry->mux[mu];
    geometry->muy[geometry->Nrays/2 + mu]  = -geometry->muy[mu];
    geometry->wmu[geometry->Nrays/2 + mu]  =  geometry->wmu[mu];
  }

  /* --- Since muz is always positive we can fill it here -- -------- */

  for (mu = 0;  mu < geometry->Nrays;  mu++) {
    geometry->muz[mu] = sqrt(1.0 - (SQ(geometry->mux[mu]) + 
				    SQ(geometry->muy[mu])));
  }
}
/* ------- end ---------------------------- getAngleQuadr.c---------- */
