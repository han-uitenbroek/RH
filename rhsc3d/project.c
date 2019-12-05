/* ------- file: -------------------------- vproject.c --------------

       Version:       rh2.0, 3-D Cartesian, short characteristics
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Wed Nov 17 08:58:43 2010 --

       --------------------------                      ----------RH-- */

/* --- Calculate projected velocity along ray mu at depth k 

   See: A. van Ballegooijen in "Numerical Radiative Transfer",
        ed. W. Kalkofen. p282., eq. 2.1
	                                --             -------------- */

#include <math.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "geometry.h"


/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern Geometry geometry;


/* ------- begin -------------------------- vproject.c -------------- */

double vproject(int k, int mu)
{
  return geometry.mux[mu]*geometry.vx[k] +
         geometry.muy[mu]*geometry.vy[k] + 
         geometry.muz[mu]*geometry.vz[k];
}
/* ------- end ---------------------------- vproject.c -------------- */

/* ------- begin -------------------------- Bproject.c -------------- */

void Bproject()
{
  /* --- Compute the cosine of gamma, the angle of B with the line of
         sight, and the sine and cosine of 2*chi, the angle of B
         with e1 --                                    -------------- */

  register int k, mu;

  double bx, by, bz, b1, b2, b3, csc_theta, sin_gamma;

  atmos.cos_gamma = matrix_double(atmos.Nrays, atmos.Nspace);
  atmos.cos_2chi  = matrix_double(atmos.Nrays, atmos.Nspace);
  atmos.sin_2chi  = matrix_double(atmos.Nrays, atmos.Nspace);

  for (mu = 0;  mu < atmos.Nrays;  mu++) {
    if (geometry.muz[mu] == 1.0) {
      for (k = 0;  k < atmos.Nspace;  k++) {
	atmos.cos_gamma[mu][k] = cos(atmos.gamma_B[k]);
	atmos.cos_2chi[mu][k]  = cos(2.0 * atmos.chi_B[k]);
	atmos.sin_2chi[mu][k]  = sin(2.0 * atmos.chi_B[k]);
      }
    } else {
      csc_theta = 1.0 / sqrt(1.0 - SQ(geometry.muz[mu]));
      for (k = 0;  k < atmos.Nspace;  k++) {
	sin_gamma = sin(atmos.gamma_B[k]);
	bx = sin_gamma * cos(atmos.chi_B[k]);
	by = sin_gamma * sin(atmos.chi_B[k]);
	bz = cos(atmos.gamma_B[k]);
	
	b3 = geometry.mux[mu]*bx + geometry.muy[mu]*by +
	  geometry.muz[mu]*bz;
	b1 = csc_theta * (bz - geometry.muz[mu]*b3);
	b2 = csc_theta * (geometry.muy[mu]*bx - geometry.mux[mu]*by);

	atmos.cos_gamma[mu][k] = b3;
	atmos.cos_2chi[mu][k]  = (SQ(b1) - SQ(b2)) / (1.0 - SQ(b3));
	atmos.sin_2chi[mu][k]  = 2.0 * b1*b2 / (1.0 - SQ(b3));
      }
    }
  }
}
/* ------- end ---------------------------- Bproject.c -------------- */
