/* ------- file: -------------------------- vproject.c --------------

       Version:       rh2.0, 1-D spherically symmetric
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Wed Feb 12 13:36:22 2003 --

       --------------------------                      ----------RH-- */

/* --- Calculate projected velocity along ray mu at depth k -- ------ */

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "geometry.h"


/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

extern Geometry geometry;


/* ------- begin -------------------------- vproject.c -------------- */

double vproject(int k, int mu)
{
  return geometry.rays[mu].xmu[k] * geometry.vel[k];
}
/* ------- end ---------------------------- vproject.c -------------- */

/* ------- begin -------------------------- Bproject.c -------------- */

void Bproject()
{
  /* --- Since we don't do magnetic fields in spherical geometry
         this routine does nothing --                  -------------- */

  return;
}
/* ------- end ---------------------------- Bproject.c -------------- */
