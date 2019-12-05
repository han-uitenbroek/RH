/* ------- file: -------------------------- riiplane.c --------------

       Version:       rh2.0, 1-D plane-parallel
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Thu Jan 22 11:02:49 2004 --

       --------------------------                      ----------RH-- */

/* --- Evaluates angle-dependent general redistribution function RII for
       ordinary and cross redistribution of lines with sharp lower and
       broadened upperlevels.

  See: - D. G. Hummer 1962, MNRAS, 125, 21-37.

       - I. Huben\'y, 1982, JQSRT 27, 593.


       The azimuthally averaged angle-dependent redistribution
       function is used in plane-parallel atmospheres.

  See: - R.W. Milkey, R.A. Shine, and D. Mihalas 1975, ApJ 202, 250-258.

       - J.H.H.J. Bruls and S.K. Solanki 1997, A&A 325, 1179-1198.


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


#define N_GAUSS_QUADR  8


/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

extern Geometry geometry;


/* ------- begin -------------------------- RII.c ------------------- */

double RII(double v_emit, double v_abs, double adamp, int mu1, int mu2)
{
  const char routineName[] = "RII";
  register int n;
  static bool_t initialize = TRUE;
  static double xg0[N_GAUSS_QUADR], wg0[N_GAUSS_QUADR];

  double theta_min, theta_plus, vmin, vplus, theta1, theta2,
    xg[N_GAUSS_QUADR], wg[N_GAUSS_QUADR], rii, mu12, muz1, muz2,
         sin_theta, cos_theta, sec_theta, csc_theta2;

  muz1 = geometry.muz[mu1];
  muz2 = geometry.muz[mu2];
  mu12 = 1.0 - SQ(muz1) - SQ(muz2);

  theta1     = acos(muz1);
  theta2     = acos(muz2);
  theta_min  = fabs(theta1 - theta2);
  theta_plus = theta1 + theta2; 

  vmin  = (v_emit - v_abs) / 2.0;
  vplus = (v_emit + v_abs) / 2.0;

  /* --- Get Gauss-Legendre integration weights --     ------------- */

  if (initialize) {
    GaussLeg(0.0, 1.0, xg0, wg0, N_GAUSS_QUADR);
    initialize = FALSE;
  }
  /* --- Integration over azimuthal angle parametrized as theta - -- */

  rii = 0.0;
  for (n = 0;  n < N_GAUSS_QUADR;  n++) {
    xg[n] = theta_min + (theta_plus - theta_min) * xg0[n];
    wg[n] = (theta_plus - theta_min) * wg0[n];

    sin_theta  = sin(xg[n]);
    cos_theta  = cos(xg[n]);
    sec_theta  = sqrt(2.0 / (1.0 + cos_theta));
    csc_theta2 = 2.0 / (1.0 - cos_theta);

    rii += (1.0 + SQ(cos_theta)) * exp(-SQ(vmin) * csc_theta2) *
      Voigt(adamp*sec_theta, vplus*sec_theta, NULL, HUI_ETAL) /
      sqrt(mu12 - SQ(cos_theta) + 2.0*muz1*muz2*cos_theta) * wg[n];
  }

  return 3.0 / (16.0 * SQ(PI)) * rii;
}
/* ------- end ---------------------------- RII.c ------------------- */
