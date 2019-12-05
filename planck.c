/* ------- file: -------------------------- planck.c ----------------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Wed Jan  4 17:46:57 2017 --

       --------------------------                      ----------RH-- */

/* --- Evaluate the Planck function at wavelength lambda.

       Input:  Nspace         -- dimension of temperature array.
               T[Nspace]      -- array of temperatures [K].
               lambda         -- wavelength [nm].

       Output: Bnu[Nspace]    -- array of Planck function values
                                 [J m^-2 s^-1 Hz^-1 sr^-1].
 *     --                                              -------------- */

 
#include <math.h>

#include "rh.h"
#include "atom.h"
#include "spectrum.h"
#include "constant.h"


#define MAX_EXPONENT  150.0

/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */


/* ------- begin -------------------------- Planck.c ---------------- */

void Planck(int Nspace, double *T, double lambda, double *Bnu)
{
  register int k;

  double hc_kla, twohnu3_c2, hc_Tkla;

  hc_kla     = (HPLANCK * CLIGHT) / (KBOLTZMANN * NM_TO_M * lambda);
  twohnu3_c2 = (2.0*HPLANCK*CLIGHT) / CUBE(NM_TO_M * lambda);

  for (k = 0;  k < Nspace;  k++) {
    hc_Tkla = hc_kla/T[k];
    if (hc_Tkla <= MAX_EXPONENT)
      Bnu[k] = twohnu3_c2 / (exp(hc_Tkla) - 1.0);
    else
      Bnu[k] = 0.0;
  }
}
/* ------- end ---------------------------- Planck.c ---------------- */
