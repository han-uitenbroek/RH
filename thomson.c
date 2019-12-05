/* ------- file: -------------------------- thomson.c ---------------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Fri Feb  9 12:55:36 2001 --

       --------------------------                      ----------RH-- */

/* --- Thomson scattering off free electrons (non-relativistic).

  See: Mihalas (1978), p. 106

       Wavelengths are given in nm, densities in m^-3, opacities in m^2,
       and emissivities in J s^-1 Hz^-1 sr^-1.
       --                                              -------------- */
 
#include <math.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "constant.h"
#include "background.h"

/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

extern Atmosphere atmos;


/* ------- begin -------------------------- Thomson.c --------------- */

void Thomson(double *chi)
{
  register int k;

  double sigma = 8.0*PI/3.0 * pow(Q_ELECTRON/(sqrt(4.0*PI*EPSILON_0) *
					      (sqrt(M_ELECTRON)*CLIGHT)), 4);

  for (k = 0;  k < atmos.Nspace;  k++)
    chi[k] = atmos.ne[k] * sigma;
}
/* ------- end ---------------------------- Thomson.c --------------- */
