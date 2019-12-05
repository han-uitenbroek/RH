/* ------- file: -------------------------- geometry.c --------------

       Version:       rh2.0, 2-D Cartesian
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Tue Oct 16 11:30:43 2001 --

       --------------------------                      ----------RH-- */

/* --- Store the up- and downwind weights and indices for interpolation
       in Geometry structure (defined in geometry.h).
       Version for rectangular 2-D grid.

       The Geometry:


                                               || || Itop
           ^ z                                 \/ \/
           |                      
       TOP ========================================================= k=0
           =   |       |       |    ^ up |              |          =   
           =   |       |       |   /     |              |          =
           =----------------------.--------------------------------=
           =   |       |       | / k-1   |              |          =  
           =   |       |     k |/        |              |          =  
           =-------------------.-----------------------------------=
           =   |       |l-1   /|l        | l+1          |          =  
      LEFT =   |       |     / |         |              |          = RIGHT
           =   |       |    /  |         |              |dz        =    
           =   |       |   /   |         |              |    dx    =    
           =--------------.----------------------------------------=
           =   |       | / k+1 |         |              |          =  
           =   |       |/down  |         |              |          =  
    BOTTOM =========================================================-> x
                                    /\ /\                               
                                    || ||  Ibottom             


       Conventions:

         -- mu[xyz] is the cosine of the angle of a ray with the
	    positive [xyz]-axis. (Note: mux^2 + muy^2 + muz^2 = 1.0)

       --                                              -------------- */


#include <math.h>
#include <stdlib.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "geometry.h"
#include "inputs.h"
#include "statistics.h"
#include "error.h"


/* --- Function prototypes --                          -------------- */

void setStencil(Geometry *geometry, Stencil *stencil, 
		enum direction direction, enum intersect intersect,
		double ds, double fraction, int l, int k,
		enum triplet triplet, enum SC_order order);

void setVertStencil(Geometry *geometry, Stencil *stencil,
		    enum direction direction, double ds, int l, int k);

void fillLongChar(Geometry *geometry, int mu, int k, enum sweep sweep);


/* --- Global variables --                             -------------- */

extern char messageStr[];
extern InputData input;

/* ------- begin -------------------------- fillMesh.c -------------- */

void fillMesh(Geometry *geometry)
{
  register int mu, k, l;

  enum    triplet vtriplet, htriplet;
  enum    SC_order horder, vorder;
  enum    boundary hboundary = geometry->hboundary; 
  int     Nx = geometry->Nx, Nz = geometry->Nz, Nrays = geometry->Nrays,
          count, Nlongchar = 0, Nstart, Nend;
  double  mux, muz, dx, dz, fraction;
  Stencil  *st;
  LongChar *lc;

  if (input.Eddington) return;

  getCPU(2, TIME_START, NULL);

  /* --- Allocate space for interpolation stencils --  -------------- */

  geometry->stencil = (Stencil **) malloc(Nrays * sizeof(Stencil *));
  for (mu = 0;  mu < Nrays;  mu++) {
    geometry->stencil[mu] = (Stencil *) malloc(Nx*Nz * sizeof(Stencil));
  }

  /* --- If the horizontal boundary conditions are PERIODIC allocate
         space for possible long characteristics --    -------------- */

  if (hboundary == PERIODIC) {
    geometry->longchar = (LongChar **) malloc(Nrays * sizeof(LongChar *));
    for (mu = 0;  mu < Nrays;  mu++)
      geometry->longchar[mu] = (LongChar *) malloc(2*Nz * sizeof(LongChar));
  } else
    geometry->longchar = NULL;

  horder = SC_QUADRATIC;

  /* --- For each direction ray go through the x-z grid and store the
         stencil of interpolation indices and coefficients in the
         Geometry structure --                         -------------- */

  for (mu = 0;  mu < Nrays;  mu++) {
    mux = geometry->mux[mu];
    muz = geometry->muz[mu];

    /* --- Distinguish between mux > 0.0, mu = 0.0, and mux < 0.0.
           Rays with mux > 0.0 run from right to left in the downsweep,
           ie. with decreasing l. --                   -------------- */

    if (mux > 0.0) {

      /* --- First evaluate the UPWIND interpolation weigths and
             indices for the DOWN sweep. These will be the DOWNWIND
             ones in the reverse (UP) sweep of the formal solution
             --                                        -------------- */

      vtriplet = LOWER_TRIPLET;
      for (k = 1;  k < Nz;  k++) {
	dz = geometry->dz[k-1];
        if (k == 1) vorder = SC_LINEAR;  else  vorder = SC_QUADRATIC;

	count  = 0;
	if (hboundary == PERIODIC) Nend = Nx;  else  Nend = Nx - 1;
	for (l = 0;  l < Nend;  l++) {
	  st = &geometry->stencil[mu][k*Nx + l];

	  dx = geometry->dx[l];
	  if (hboundary == FIXED  &&  l == Nx-2)
            htriplet = LOWER_TRIPLET;
	  else
	    htriplet = UPPER_TRIPLET;
	  fraction = mux*dz / (muz*dx);

	  /* --- Characteristic intersects horizontal grid line? -- - */

	  if (fraction <= 1.0) {
	    count++;
            setStencil(geometry, st, UPWIND, HORIZONTAL, dz/muz, 
		       fraction, l, k-1, htriplet, horder);
	  } else {
	    setStencil(geometry, st, UPWIND, VERTICAL, dx/mux,
		       1.0 - 1.0/fraction, l+1, k-1, vtriplet, vorder);
	  }
	}
	/* --- If no horizontal intersection occurred use long
	       characteristic, otherwise find proper starting point - */
	
	if (geometry->longchar != NULL) {
	  lc = &geometry->longchar[mu][DOWN*Nz + k];
	  if (count == 0) {
	    lc->lstart = Nx-1;
	    fillLongChar(geometry, mu, k, DOWN);
	    Nlongchar++;
	  } else {
	    lc->Nlc    = 0;
	    lc->lstart = 0;
	    st = &geometry->stencil[mu][k*Nx + lc->lstart];
	    while(st->intersect[UPWIND] == VERTICAL) {
	      lc->lstart++;
	      st++;
	    }
	  }
	}
      }
      /* --- Evaluate the DOWNWIND interpolation weights and indices
             for the DOWN sweep. These will be the UPWIND ones in
             the reverse (UP) sweep. --                  ------------ */

      vtriplet = UPPER_TRIPLET;
      for (k = 0;  k < Nz-1;  k++) {
        dz = geometry->dz[k];
        if (k == Nz - 2) vorder = SC_LINEAR;  else  vorder = SC_QUADRATIC;

	count = 0;
	if (hboundary == PERIODIC) Nstart = 0;  else  Nstart = 1;
	for (l = Nstart;  l < Nx;  l++) {
	  st = &geometry->stencil[mu][k*Nx + l];

	  if (hboundary == FIXED) {
	    dx = geometry->dx[l-1];
	    if (l == 1)
	      htriplet = UPPER_TRIPLET;
	    else
	      htriplet = LOWER_TRIPLET;
	  } else {
	    dx = geometry->dx[MODULO(l-1, Nx)];
            htriplet = LOWER_TRIPLET;
	  }
	  fraction = mux*dz / (muz*dx);

	  /* --- Characteristic intersects horizontal grid line? -- - */

	  if (fraction <= 1.0) {
	    count++;
	    setStencil(geometry, st, DOWNWIND, HORIZONTAL, dz/muz,
		       1.0 - fraction, l-1, k+1, htriplet, horder);
	  } else {
	    setStencil(geometry, st, DOWNWIND, VERTICAL, dx/mux,
		       1.0/fraction, l-1, k, vtriplet, vorder);
	  }
	}
	/* --- If no horizontal intersection occurred use long
	       characteristic, otherwise find proper starting point - */

	if (geometry->longchar != NULL) {
	  lc = &geometry->longchar[mu][UP*Nz + k];
	  if (count == 0) {
	    lc->lstart = 0;
	    fillLongChar(geometry, mu, k, UP);
	    Nlongchar++;
	  } else {
	    lc->Nlc    = 0;
	    lc->lstart = Nx - 1;
	    st = &geometry->stencil[mu][k*Nx + lc->lstart];
	    while(st->intersect[DOWNWIND] == VERTICAL) {
	      lc->lstart--;
	      st--;
	    }
	  }
	}
      }
    } else if (mux == 0.0) {

      /* --- Case for purely vertical rays --          -------------- */

      for (k = 1;  k < Nz;  k++) {
	dz = geometry->dz[k-1];
        for (l = 0;  l < Nx;  l++) {
	  st = &geometry->stencil[mu][k*Nx + l];
	  setVertStencil(geometry, st, UPWIND, dz, l, k-1);
	}
	if (geometry->hboundary == PERIODIC) {
	  lc = &geometry->longchar[mu][DOWN*Nz + k];
	  lc->Nlc    = 0;
	  lc->lstart = 0;
	}
      }
      for (k = 0;  k < Nz-1;  k++) {
	dz = geometry->dz[k];
        for (l = 0;  l < Nx;  l++) {
	  st = &geometry->stencil[mu][k*Nx + l];
	  setVertStencil(geometry, st, DOWNWIND, dz, l, k+1);
	}
	if (geometry->hboundary == PERIODIC) {
	  lc = &geometry->longchar[mu][UP*Nz + k];
	  lc->Nlc    = 0;
	  lc->lstart = 0;
	}
      }
    } else {

      /* --- Rays with mux < 0.0 run from left to right in the
             downsweep, ie. with increasing l.
             First evaluate the UPWIND interpolation weigths and
             indices for the DOWN sweep. These will be the DOWNWIND
             ones in the reverse (UP) sweep of the formal solution
             --                                        -------------- */

      vtriplet = LOWER_TRIPLET;
      for (k = 1;  k < Nz;  k++) {
	dz = geometry->dz[k-1];
	if (k == 1) vorder = SC_LINEAR;  else  vorder = SC_QUADRATIC;

	count  = 0;
	if (hboundary == PERIODIC) Nstart = 0;  else  Nstart = 1;
	for (l = Nstart;  l < Nx;  l++) {
	  st = &geometry->stencil[mu][k*Nx + l];

	  if (hboundary == FIXED) {
	    dx = geometry->dx[l-1];
	    if (l == 1)
	      htriplet = UPPER_TRIPLET;
            else
	      htriplet = LOWER_TRIPLET;
	  } else {
	    dx = geometry->dx[MODULO(l-1, Nx)];
	    htriplet = LOWER_TRIPLET;
	  }
	  fraction = -mux*dz / (muz*dx);

	  /* --- Characteristic intersects horizontal grid line? -- - */

	  if (fraction <= 1.0) {
	    count++;
	    setStencil(geometry, st, UPWIND, HORIZONTAL, dz/muz,
		       1.0 - fraction, l-1, k-1, htriplet, horder);
	  } else {
	    setStencil(geometry, st, UPWIND, VERTICAL, -dx/mux,
		       1.0 - 1.0/fraction, l-1, k-1, vtriplet, vorder);
	  }
	}
	/* --- If no horizontal intersection occurred use long
	       characteristic, otherwise find proper starting point - */
	
	if (geometry->longchar != NULL) {
	  lc = &geometry->longchar[mu][DOWN*Nz + k];
	  if (count == 0) {
            lc->lstart = 0;
	    fillLongChar(geometry, mu, k, DOWN);
	    Nlongchar++;
	  } else {
	    lc->Nlc    = 0;
	    lc->lstart = Nx-1;
	    st = &geometry->stencil[mu][k*Nx + lc->lstart];
	    while(st->intersect[UPWIND] == VERTICAL) {
	      lc->lstart--;
	      st--;
	    }
	  }
	}
      }
      /* --- Evaluate the DOWNWIND interpolation weights and indices
             for the DOWN sweep. These will be the UPWIND ones in
             the reverse (UP) sweep. --                  ------------ */

      vtriplet = UPPER_TRIPLET;
      for (k = 0;  k < Nz-1;  k++) {
        dz = geometry->dz[k];
	if (k == Nz - 2) vorder = SC_LINEAR;  else  vorder = SC_QUADRATIC;

	count  = 0;
	if (hboundary == PERIODIC) Nend = Nx;  else  Nend = Nx - 1;
	for (l = 0;  l < Nend;  l++) {
	  st = &geometry->stencil[mu][k*Nx + l];

	  dx = geometry->dx[l];
	  if (hboundary == FIXED  &&  l == Nx - 2)
	    htriplet = LOWER_TRIPLET;
	  else
            htriplet = UPPER_TRIPLET;
	  fraction = -mux*dz / (muz*dx);

	  /* --- Characteristic intersects horizontal grid line? -- - */

	  if (fraction <= 1.0) {
	    count++;
	    setStencil(geometry, st, DOWNWIND, HORIZONTAL, dz/muz,
		       fraction, l, k+1, htriplet, horder);
	  } else {
	    setStencil(geometry, st, DOWNWIND, VERTICAL, -dx/mux,
		       1.0/fraction, l+1, k, vtriplet, vorder);
	  }
	}
	/* --- If no horizontal intersection occurred use long
	       characteristic, otherwise find proper starting point - */

	if (geometry->longchar != NULL) {
	  lc =&geometry->longchar[mu][UP*Nz + k];
	  if (count == 0) {
	    lc->lstart = Nx-1;
	    fillLongChar(geometry, mu, k, UP);
	    Nlongchar++;
	  } else {
	    lc->Nlc    = 0;
	    lc->lstart = 0;
	    st = &geometry->stencil[mu][k*Nx + lc->lstart];
	    while(st->intersect[DOWNWIND] == VERTICAL) {
	      lc->lstart++;
	      st++;
	    }
	  }
	}
      }
    }
  }

  sprintf(messageStr, "\n-Used %d Long Characteristics\n\n", Nlongchar);
  Error(MESSAGE, NULL, messageStr);
  getCPU(2, TIME_POLL, "Geometry");
}
/* ------- end ---------------------------- fillMesh.c -------------- */

/* ------- begin -------------------------- setVertStencil.c -------- */

void setVertStencil(Geometry *geometry, Stencil *stencil,
		    enum direction direction, double ds, int l, int k)
{
  const char routineName[] = "setVertStencil";
  register int n;

  int Nx = geometry->Nx;

  stencil->intersect[direction] = HORIZONTAL;
  stencil->ds[direction]        = ds;
  stencil->fraction[direction]  = 0.0;
  stencil->order[direction]     = SC_LINEAR;
  stencil->triplet[direction]   = LOWER_TRIPLET;

  for (n = 0;  n < 3;  n++)
    stencil->index[direction][n] = k*Nx + l;

  stencil->coeff[direction][0] = 1.0;
  stencil->coeff[direction][1] = 0.0;
  stencil->coeff[direction][2] = 0.0;
}
/* ------- end ---------------------------- setVertStencil.c -------- */

/* ------- begin -------------------------- setStencil.c ------------ */

void setStencil(Geometry *geometry, Stencil *stencil, 
		enum direction direction, enum intersect intersect,
		double ds, double fraction, int l, int k, 
		enum triplet triplet, enum SC_order order)
{
  const char routineName[] = "setStencil";
  register int n;

  int    Nx = geometry->Nx, index[3];
  double d0, d1, d2, ds0, ds1, ds2;

  if (fraction < 0.0) {
    Error(ERROR_LEVEL_2, routineName,
	  "Fractional index for interpolation should be > 0.0\n");
  }
  switch (direction) {
  case UPWIND:   break;
  case DOWNWIND: break;
  default:
    Error(ERROR_LEVEL_2, routineName,
	  "Invalid direction index for Stencil");
  }
  stencil->intersect[direction] = intersect;
  stencil->ds[direction]        = ds;
  stencil->fraction[direction]  = fraction;
  stencil->triplet[direction]   = triplet;
  stencil->order[direction]     = order;

  if (order == SC_LINEAR) {
    if (intersect == HORIZONTAL) {
      if (geometry->hboundary == PERIODIC) {
	stencil->index[direction][0] = k*Nx + MODULO(l, Nx);
	stencil->index[direction][1] = k*Nx + MODULO(l+1, Nx);
      } else {
	stencil->index[direction][0] = k*Nx + l;
	stencil->index[direction][1] = k*Nx + l+1;
      }
    } else {
      if (geometry->hboundary == PERIODIC) l = MODULO(l, Nx);
      stencil->index[direction][0] = k*Nx + l;
      stencil->index[direction][1] = (k+1)*Nx + l;
    }

    stencil->coeff[direction][0] = (1.0 - fraction);
    stencil->coeff[direction][1] = fraction;
  } else {

    if (triplet == LOWER_TRIPLET) {
      if (intersect == HORIZONTAL) l--;  else  k--;
    }
    if (intersect == HORIZONTAL) {
      for (n = 0;  n < 3;  n++) {
	if (geometry->hboundary == PERIODIC)
	  index[n] = MODULO(l + n, Nx);
	else
          index[n] = l + n;
	stencil->index[direction][n] = k*Nx + index[n];
      }
      ds0 = geometry->dx[index[0]];
      ds1 = geometry->dx[index[1]];
    } else {
      if (geometry->hboundary == PERIODIC) l = MODULO(l, Nx);
      for (n = 0;  n < 3;  n++) {
        index[n] = k + n;
	stencil->index[direction][n] = index[n]*Nx + l;
      }
      ds0 = geometry->dz[index[0]];
      ds1 = geometry->dz[index[1]];
    }
    ds2 = ds0 + ds1;

    if (triplet == UPPER_TRIPLET) {
      d0 = fraction * ds0;
      d1 = ds0 - d0;
    }
    else {
      d1 = -fraction * ds1;
      d0 =  ds0 - d1;
    }
    d2 = d1 + ds1;

    stencil->coeff[direction][0] =  d1*d2 / (ds0*ds2);
    stencil->coeff[direction][1] =  d0*d2 / (ds0*ds1);
    stencil->coeff[direction][2] = -d0*d1 / (ds2*ds1);
  }
}
/* ------- end ---------------------------- getcoeffs.c ------------- */

/* ------- begin -------------------------- fillLongChar.c ---------- */

void fillLongChar(Geometry *geometry, int mu, int k, enum sweep sweep)
{
  register int l, n;

  enum triplet htriplet, vtriplet;
  enum SC_order vorder;
  int    Nx = geometry->Nx, Nz = geometry->Nz, dl, k_upwind, count,
         l_intersect;
  double mux, muz, dx_total, dx_intersect, dz, dz_total, ds_total,
         ds, fraction;
  LongChar *lc;

  /* --- Fills the appropriate LongChar structure of the interpolation
         mesh with indices and quadratic coefficients to solve along a
         long characteristic when, with periodic or even boundary
         conditions, no shortcharacteristic intersects a horizontal grid
         line for a given k. --                        -------------- */

  lc  = geometry->longchar[mu]+(sweep*Nz + k);
  mux = geometry->mux[mu];
  muz = geometry->muz[mu];

  l = lc->lstart;
  if (sweep == DOWN) {
    vtriplet = LOWER_TRIPLET;
    if (k == 1) vorder = SC_LINEAR;  else  vorder = SC_QUADRATIC;
    k_upwind = k-1;
    dz = geometry->dz[k_upwind];
    if (mux >= 0.0) {
      htriplet = UPPER_TRIPLET;
      dl = 1;
    } else {
      htriplet = LOWER_TRIPLET;
      dl = -1;
    }
  } else {
    vtriplet = UPPER_TRIPLET;
    if (k == Nz-2) vorder = SC_LINEAR;  else  vorder = SC_QUADRATIC;
    k_upwind = k+1;
    dz = geometry->dz[k];
    if (mux >= 0.0) {
      htriplet = LOWER_TRIPLET;
      dl = -1;
    } else {
      htriplet = UPPER_TRIPLET;
      dl = 1;
    }
  }
  /* --- Find first intersection with horizontal grid line -- ------- */

  count = 0;
  dx_total = 0.0;
  dx_intersect = dz/muz * fabs(mux);

  if (dl > 0) {
    while (dx_total < dx_intersect) {
      dx_total += geometry->dx[l];
      l = MODULO(l + 1, Nx);
      count++;
    }
    l = MODULO(l - 1, Nx);
  } else {
    while (dx_total < dx_intersect) {
      l = MODULO(l - 1, Nx);
      dx_total += geometry->dx[l];
      count++;
    }
  }
  /* --- Allocate space for long characteristic indices, interpolation
         coefficients, and pathlength.

   Note: Interpolation data is stored in an array of Stencil structures
         so that the setStencil function can be used to determine the
         interpolation coefficients, and Quadr can be used to do the 
         actual interpolation. Only the UPWIND half of each stencil is
         used because the long characteristic is solved only in one
         direction, towards the root. --               -------------- */

  lc->Nlc     = count;
  lc->stencil = (Stencil *) malloc(count * sizeof(Stencil));

  /* --- Store the interpolation in the horizontal grid line -- ----- */

  if (dl > 0) {
    fraction = 1.0 - (dx_total - dx_intersect)/geometry->dx[l];
    ds       = fraction * geometry->dx[l] / fabs(mux);
  } else {
    fraction = (dx_total - dx_intersect)/geometry->dx[l];
    ds       = (1.0 - fraction)*geometry->dx[l] / fabs(mux);
  }
  setStencil(geometry, &lc->stencil[0], UPWIND, HORIZONTAL, ds,
	     fraction, l, k_upwind, htriplet, SC_QUADRATIC);

  /* --- Then segment for segment the vertical crossings -- --------- */

  ds_total = 0.0;
  for (n = 1;  n < count;  n++) {
    ds_total += lc->stencil[n-1].ds[UPWIND];
    dz_total  = ds_total * muz;
    if (dl > 0) l_intersect = l;  else  l_intersect = MODULO(l + 1, Nx);

    l = MODULO(l - dl, Nx);
    ds = geometry->dx[l] / fabs(mux);
    if (sweep == DOWN) {
      fraction = dz_total / dz;
      setStencil(geometry, &lc->stencil[n], UPWIND, VERTICAL, ds,
		 fraction, l_intersect, k_upwind, vtriplet, vorder);
    } else {
      fraction = 1.0 - dz_total / dz;
      setStencil(geometry, &lc->stencil[n], UPWIND, VERTICAL, ds,
		 fraction, l_intersect, k, vtriplet, vorder);
    }
  }
}
/* ------- end ---------------------------- fillLongChar.c ---------- */
