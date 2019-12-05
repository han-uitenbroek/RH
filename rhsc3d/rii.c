/* ------- file: -------------------------- rii.c -------------------

       Version:       rh2.0, 3-D Cartesian
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Wed Oct 18 16:08:12 2000 --

       --------------------------                      ----------RH-- */

/* --- Evaluates angle-dependent general redistribution function RII for
       ordinary and cross redistribution of lines with sharp lower and
       broadened upperlevels.

  See: - D. G. Hummer 1962, MNRAS, 125, 21-37.

       - I. Huben\'y, 1982, JQSRT 27, 593.


       Parameters:

         v_abs  -- Frequency absorbed photon (Doppler units).
         v_emit -- Frequency emitted photon (Doppler units).
         adamp  -- Damping parameter line.
         mu1,2  -- Index of incoming and outgoing ray, respectively.

       Note: The dipole phase function is used:

         g(n, n') = 3/4 (1 + cos^2(theta)).

       --                                              -------------- */

#include <math.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "geometry.h"
#include "constant.h"
#include "error.h"


/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

extern Geometry geometry;


/* ------- begin -------------------------- RII.c ------------------- */

double RII(double v_emit, double v_abs, double adamp, int mu1, int mu2)
{
  const char routineName[] = "RII";

  double vmin, vplus, rii, mu12,
         sin_theta, cos_theta, sec_theta, csc_theta2;

  if (mu1 == mu2) return 0.0;

  /* --- cos_theta is the cosine of the acute angle between incoming
         and outgoing photons. --                      -------------- */

  cos_theta = geometry.mux[mu1] * geometry.mux[mu2] +
              geometry.muy[mu1] * geometry.muy[mu2] + 
              geometry.muz[mu1] * geometry.muz[mu2];
  if (cos_theta < 0.0) cos_theta = -cos_theta;

  sin_theta  = sqrt(1.0 - SQ(cos_theta));
  sec_theta  = sqrt(2.0 / (1.0 + cos_theta));
  csc_theta2 = 2.0 / (1.0 - cos_theta);

  vmin  = (v_emit - v_abs) / 2.0;
  vplus = (v_emit + v_abs) / 2.0;

  rii = ((1.0 + SQ(cos_theta))/sin_theta) * exp(-SQ(vmin)*csc_theta2) *
    Voigt(adamp * sec_theta, vplus * sec_theta, NULL, HUI_ETAL);

  return 3.0 / (16.0 * SQ(PI)) * rii;
}
/* ------- end ---------------------------- RII.c-------------------- */
