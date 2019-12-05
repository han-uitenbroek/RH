/* ------- file: -------------------------- geometry.c --------------

       Version:       rh2.0, 1-D spherically symmetric
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Thu Aug 16 07:44:51 2001 --

       --------------------------                      ----------RH-- */

/* --- Figure out the geometry of the different rays in spherical
       geometry.

  See: A. Nordlund (1984), in "Methods in radiative transfer",
       W. Kalkofen (ed.), p. 211-233.
       --                                              -------------- */

 
#include <stdlib.h>
#include <math.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "geometry.h"
#include "statistics.h"

/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

extern Atmosphere atmos;

/* ------- begin -------------------------- Geometry.c -------------- */

void getRays(Geometry *geometry)
{
  register int mu, k, n;

  int     Nradius = geometry->Nradius, ncore, Nray0;
  double *pcore, xmu, *r = geometry->r;
  Ray *ray;

  getCPU(2, TIME_START, NULL);

  /* --- PTOP is the depth index of the outermost shell which
         is allowed to have a grazing ray --           -------------- */

  geometry->Nrays = geometry->Nradius - PTOP + geometry->Ncore;
  atmos.Nrays     = geometry->Nrays;
  geometry->rays  = (Ray *) malloc(geometry->Nrays * sizeof(Ray));

  /* --- Figure out impact parameters for core rays, spaced equidistant
         in mu --                                      -------------- */

  pcore = (double *) malloc(geometry->Ncore * sizeof(double));
  for (n = 0;  n < geometry->Ncore;  n++) {
    xmu = (n + 1) / ((double) geometry->Ncore);
    pcore[n] = sqrt(1.0 - SQ(xmu)) * (r[Nradius - 1] + geometry->Radius);
  }
  /* --- Go through all the rays and fill ray structures -- --------- */

  ncore = 0;
  for (mu = 0;  mu < geometry->Nrays;  mu++) {
    ray = &geometry->rays[mu];
    if (mu < (geometry->Nradius - PTOP)) {
      ray->type = SHELL_RAY;
      ray->Ns   = mu + PTOP + 1;
      ray->p    = r[mu + PTOP];
    } else {
      ray->type = CORE_RAY;
      ray->Ns   = geometry->Nradius;
      ray->p    = pcore[ncore++] - geometry->Radius;
    }
    ray->s   = (double *) malloc(ray->Ns * sizeof(double));
    ray->xmu = (double *) malloc(ray->Ns * sizeof(double));
    ray->wmu = (double *) malloc(ray->Ns * sizeof(double));
    for (k = 0;  k < ray->Ns;  k++) {
      ray->s[k]   = sqrt((r[k] - ray->p) *
			 (2.0*geometry->Radius + r[k] + ray->p));
      ray->xmu[k] = ray->s[k] / (geometry->Radius + r[k]);
    }
  }
  /* --- Collect the angular integration weights --    -------------- */

  for (k = 0;  k < geometry->Nradius;  k++) {
    if (k < PTOP) {
      ray = &geometry->rays[0];
      ray->wmu[k] = 0.5*((ray+1)->xmu[k] - ray->xmu[k]) +
        ray->xmu[k];
    } else {
      ray = &geometry->rays[k - PTOP];
      ray->wmu[k] = 0.5*((ray+1)->xmu[k] - ray->xmu[k]);
    }
    ray = ray + 1;
    for (mu = MAX(0, k-PTOP)+1;  mu < geometry->Nrays-1;  mu++) {
      ray->wmu[k] = 0.5*((ray+1)->xmu[k] - (ray-1)->xmu[k]);
      ray++;
    }
    ray->wmu[k] = 0.5*(ray->xmu[k] - (ray-1)->xmu[k]);
  }
  /* --- Clean up --                                   -------------- */

  free(pcore);
  getCPU(2, TIME_POLL, "Geometry");
}
/* ------- end ---------------------------- Geometry.c -------------- */
