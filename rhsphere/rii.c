/* ------- file: -------------------------- rii.c -------------------

       Version:       rh2.0, , 1-D spherically symmetric
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Wed Oct 18 16:25:02 2000 --

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

#include "rh.h"
#include "atom.h"
#include "error.h"


/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */


/* ------- begin -------------------------- RII.c ------------------- */

double RII(double v_emit, double v_abs, double adamp, int mu1, int mu2)
{
  const char routineName[] = "RII";

  Error(ERROR_LEVEL_2, routineName,
	"Not yet implemented for SPHERICAL_SYMMETRY");

  return 0.0;
}
/* ------- end ---------------------------- RII.c-------------------- */
