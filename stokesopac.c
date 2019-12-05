/* ------- file: -------------------------- stokesopac.c ------------

       Version:       rh1.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Mon Feb 23 15:51:02 2004 --

       --------------------------                      ----------RH-- */

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "spectrum.h"
#include "inputs.h"

/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern Spectrum spectrum;
extern InputData input;
extern char messageStr[];


/* ------- begin -------------------------- StokesK.c --------------- */

void StokesK(int nspect, int k, double chi_I, double K[4][4])
{
  register int i, j;

  ActiveSet *as;

  /* --- Return the elements of the reduced 4x4 Stokes propagation
         matrix K', which is defined as:

           =           =         =              =
           K' = (chi_c*1 + chi_l*Phi) / chi_I - 1,

	   for wavelength# nspect, spatial point k, and ray mu.

    See: Rees, Murphy, & Durrant, 1989, ApJ 339, 1093-1106.

         --                                            -------------- */
  as = &spectrum.as[nspect];

  for (j = 0;  j < 4;  j++)
    for (i = 0;  i < 4;  i++) K[j][i] = 0.0;

  if (containsPolarized(as)) {
    K[0][1] = as->chi[atmos.Nspace + k];
    K[0][2] = as->chi[2*atmos.Nspace + k];
    K[0][3] = as->chi[3*atmos.Nspace + k];

    if (input.magneto_optical) {
      K[1][2] = as->chip[2*atmos.Nspace + k];
      K[1][3] = as->chip[atmos.Nspace + k];
      K[2][3] = as->chip[k];
    }
  }
  if (atmos.backgrflags[nspect].ispolarized) {
    K[0][1] += as->chi_c[atmos.Nspace + k];
    K[0][2] += as->chi_c[2*atmos.Nspace + k];
    K[0][3] += as->chi_c[3*atmos.Nspace + k];

    if (input.magneto_optical) {
      K[1][2] += as->chip_c[2*atmos.Nspace + k];
      K[1][3] += as->chip_c[atmos.Nspace + k];
      K[2][3] += as->chip_c[k];
    }
  }
  /* --- Divide by Stokes I opacity and fill lower diagonal part -- - */

  for (j = 0;  j < 3;  j++) {
    for (i = j+1;  i < 4;  i++) {
      K[j][i] /= chi_I;
      K[i][j]  = K[j][i];
    }
  }
  /* --- Anti-symmetric magneto-optical elements --    -------------- */

  if (input.magneto_optical) {
    K[1][3] *= -1.0;
    K[2][1] *= -1.0;
    K[3][2] *= -1.0;
  }
}
/* ------- end ---------------------------- StokesK.c --------------- */
