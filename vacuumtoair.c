/* ------- file: -------------------------- vacuumtoair.c -----------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Thu Sep  8 22:17:50 2005 --

       --------------------------                      ----------RH-- */

/* --- Convert wavelengths at and above VACUUM_TO_AIR_LIMIT
       (see spectrum.h) to air wavelengths.
       An approximation to the 4th power of inverse wavenumber is used

  See: IUE Image Processing Manual Page 6-15.
       --                                              -------------- */
 
#include "rh.h"
#include "atom.h"
#include "spectrum.h"


/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */


/* ------- begin -------------------------- vacuum_to_air.c --------- */

void vacuum_to_air(int Nlambda, double *lambda_vac, double *lambda_air)
{
  register int la;

  double sqwave, reduction;

  /* --- Wavelengths should be in nm. --               -------------- */

  for (la = 0;  la < Nlambda;  la++) {
    if (lambda_vac[la] >= VACUUM_TO_AIR_LIMIT) {
      sqwave = 1.0 / SQ(lambda_vac[la]);
      reduction = 1.0 + 2.735182E-4 +
	(1.314182 + 2.76249E+4 * sqwave) * sqwave;
      lambda_air[la] = lambda_vac[la] / reduction;
    } else
      lambda_air[la] = lambda_vac[la];
  }
}
/* ------- end ---------------------------- vacuum_to_air.c --------- */

/* ------- begin -------------------------- air_to_vacuum.c --------- */

void air_to_vacuum(int Nlambda, double *lambda_air, double *lambda_vac)
{
  register int la;

  double sqwave, increase;

  /* --- Wavelengths should be in nm. --               -------------- */

  for (la = 0;  la < Nlambda;  la++) {
    if (lambda_air[la] >= AIR_TO_VACUUM_LIMIT) {
      sqwave = SQ(1.0E+07 / lambda_air[la]);
      increase = 1.0000834213E+00 +
	2.406030E+06/(1.30E+10 - sqwave) + 1.5997E+04/(3.89E+09 - sqwave);
      lambda_vac[la] = lambda_air[la] * increase;
    } else
      lambda_vac[la] = lambda_air[la];
  }
}
/* ------- end ---------------------------- air_to_vacuum.c --------- */
