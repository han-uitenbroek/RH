/* ------- file: -------------------------- anglequad.c -------------

       Version:       rh2.0, 1-D plane-parallel
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Fri Jun 17 16:18:00 2022 --

       By default the quadrature points are the zeroes of the Gauss-
       Legendre polynomials.

       There is an option for Gauss-Lobatto
       quadrature, which includes the vertical direction
         See: Abramowitz & Stegun 1972, p. 888

       --------------------------                      ----------RH-- */

/* --- Get angle quadrature. --                        -------------- */

 
#include <stdlib.h>
#include <math.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "geometry.h"


/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

extern Atmosphere atmos;


/* ------- begin -------------------------- getAngleQuad.c ---------- */

void getAngleQuad(Geometry *geometry)
{
  register int mu;

  /* --- Points and Weights for the Gauss-Lobatto quadrature -- ----- */
  
  double muz_GLOB_4[NRO_GLOB_4] = {sqrt(1.0/5.0), 1.0}; 
  double wmu_GLOB_4[NRO_GLOB_4] = {5.0/6.0, 1.0/6.0};

  double muz_GLOB_6[NRO_GLOB_6] = {sqrt(1.0/3.0 - (2.0*sqrt(7.0)/21.0)),
				   sqrt(1.0/3.0 + (2.0*sqrt(7.0)/21.0)),
				   1.0};
  double wmu_GLOB_6[NRO_GLOB_6] = {(14.0 + sqrt(7.0))/30.0,
				   (14.0 - sqrt(7.0))/30.0,
				   1.0/15.0};
  
  
  /* --- Copy the number of rays to the geometry structure.
   Note: Nrays is read into the atmos structure in readinput.c -- --- */

  geometry->Nrays = atmos.Nrays;

  geometry->mux = (double *) malloc(geometry->Nrays * sizeof(double));
  geometry->muy = (double *) malloc(geometry->Nrays * sizeof(double));
  geometry->muz = (double *) malloc(geometry->Nrays * sizeof(double));
  geometry->wmu = (double *) malloc(geometry->Nrays * sizeof(double));

  switch (atmos.angleSet.set) {
  case SET_VERTICAL:
    geometry->muz[0] = 1.0;
    geometry->wmu[0] = 1.0;
    break;
  case SET_EDDINGTON:
    geometry->muz[0] = 1.0 / sqrt(3.0);
    geometry->wmu[0] = 1.0;
    break;
  case SET_GAUSS_LOBATTO:
    switch (geometry->Nrays) {
    case NRO_GLOB_4:
      for (mu = 0;  mu < NRO_GLOB_4;  mu++) {
	geometry->muz[mu] = muz_GLOB_4[mu];
	geometry->wmu[mu] = wmu_GLOB_4[mu];
      } break;
    case NRO_GLOB_6:
       for (mu = 0;  mu < NRO_GLOB_6;  mu++) {
	geometry->muz[mu] = muz_GLOB_6[mu];
	geometry->wmu[mu] = wmu_GLOB_6[mu];
      } break;
    } break;
  default:
    GaussLeg(0.0, 1.0, geometry->muz, geometry->wmu, geometry->Nrays);
  }

  /* --- In the 1-D plane case assume that all the rays lie in the z-x
         plane and are pointed in the direction of positive x. -- --- */

  for (mu = 0;  mu < geometry->Nrays;  mu++) {
    geometry->muy[mu] = 0.0;
    geometry->mux[mu] = sqrt(1.0 - SQ(geometry->muz[mu]));
  }
}
/* ------- end ---------------------------- getAngleQuad.c ---------- */
