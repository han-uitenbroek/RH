/* ------- file: -------------------------- geometry.c --------------
 
       Version:       rh2.0, 3-D Cartesian, short characteristics
       Authors:       Han Uitenbroek   (huitenbroek@nso.edu)
                      Jorrit Leenaarts (j.leenaarts@astro.uu.nl)

       Last modified: Wed Apr 22 09:08:00 2009 --
 
       --------------------------                      ----------RH-- */
 
/* --- Store the up- and downwind interpolation bases, stencil and ds's
       in the Geometry structure (defined in geometry.h) pointed to by 
       geometry.

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
           =   |       |     / |         |              |          =
           =   |       |    /  |         |              |dz        =
           =   |       |   /   |         |              |    dx    =/dy
           =--------------.----------------------------------------=
           =   |       | / k+1 |         |              |          = y
           =   |       |/down  |         |              |          =/
    BOTTOM =========================================================-> x
                                    /\ /\
                                    || ||  Ibottom

       Currently allow only PERIODIC boundary conditions in x and y,
       and equidistant grids in the x- and y-directions. The x- and
       y grid may have different periods and spacing. The z-grid may
       be irregular.

  See: P. B. Kunasz and  L. H. Auer 1988, JQSRT 39, 67-79

       L. H. Auer, P. Fabiani Bendicho and J. Trujillo Bueno 1994,
       A&A, 292, 599-615

       At the top and bottom only a fixed boundary is allowed.

 Note: Rays with to_obs == TRUE go from BOTTOM to TOP.


       Cubic convolution can be used in the horizontal plane if keyword
       input.interpolate_3D == BICUBIC_3D. Otherwise, linear
       interpolation is used as it is in the vertical planes.

  See: R.G. Keys, 1981, in IEEE Trans. Acoustics, Speech,
       and Signal Processing, Vol. 29, pp. 1153-1160.

       --                                              -------------- */

#include <math.h>
#include <stdlib.h>

#include "rh.h"
#include "error.h"
#include "atom.h"
#include "atmos.h"
#include "geometry.h"
#include "inputs.h"
#include "statistics.h"
#include "error.h"

#define MIN_FRAC 0.1


/* --- Function prototypes --                          -------------- */

double  cubicconvol(double *f, int Nx, int *x_base, int *y_base,
		    double *x_kernel, double *y_kernel, bool_t monotonic);
int     sascend(const void *vis1, const void *vis2);
double  round(double x);
double  calcDistance(double frac1, double frac2, enum intersect plane,
		     double dx, double dy, double dz);
void    calcFrac(double *frac1, double *frac2,
		 double rx, double ry, double rz, enum intersect plane,
		 enum direction wind);
enum    intersect getIntersect(double rx, double ry, double rz);
int     getQuadrant(double mux, double muy, enum direction wind);
void    setStencil(Geometry *geometry, int mu, int k,
		   enum direction wind);


/* --- Global variables --                             -------------- */ 

extern InputData input;
extern char messageStr[];
  

/* ------- begin -------------------------- fillMesh.c -------------- */ 

void fillMesh(Geometry *geometry)
{
  const char routineName[] = "fillMesh";
  register int mu, k, l, m;
 
  int    quadrant;
  double mux, muy, muz, dz, frac1, frac2, longfracx, longfracy;
  enum intersect plane;

  getCPU(2, TIME_START, NULL);
 
  if (geometry->x_boundary != PERIODIC ||
      geometry->y_boundary != PERIODIC) {
    Error(ERROR_LEVEL_2, routineName, 
	  "Horizontal boundary conditions in x AND y must be PERIODIC");
  }
  /* --- Allocate space for interpolation stencils and longchar
         structures --                                 -------------- */

  geometry->stencil_uw =
    (Stencil **) malloc(geometry->Nrays * sizeof(Stencil *));
  geometry->stencil_dw =
    (Stencil **) malloc(geometry->Nrays * sizeof(Stencil *));

  for (mu = 0;  mu < geometry->Nrays;  mu++) {
    geometry->stencil_uw[mu] =
      (Stencil *) malloc(geometry->Nz * sizeof(Stencil));
    geometry->stencil_dw[mu] =
      (Stencil *) malloc(geometry->Nz * sizeof(Stencil));
  }
  geometry->Nlongchar = 0;
   
  /* --- Evaluate the interpolation bases and kernels for linear
         interpolation or cubic convolution. In general we sweep DOWN
         and then UP. The convention here will be that we store the
         UPWIND base for the DOWN sweep in the uw_stencil array, and the
         DOWNWIND in the dw_stencil array. For the UP sweep UPWIND and
         DOWNWIND are reversed. --                     -------------- */

  for (mu = 0;  mu < geometry->Nrays;  mu++) {
    for (k = 1;  k < geometry->Nz;  k++) {

      /* --- First, the UPWIND bases for the DOWN sweep. These will be
             the DOWNWIND bases for the final UP sweep. -- ---------- */ 

      setStencil(geometry, mu, k, UPWIND);
    }
    for (k = 0;  k < geometry->Nz - 1;  k++) {
      setStencil(geometry, mu, k, DOWNWIND);
    }
  }
  if (geometry->Nlongchar) {
    sprintf(messageStr, "\n-Used %d Long Characteristics\n\n",
	    geometry->Nlongchar);
    Error(MESSAGE, NULL, messageStr);
  }
  getCPU(2, TIME_POLL, "Geometry");
}
/* ------- end ---------------------------- fillMesh.c -------------- */

/* ------- begin -------------------------- getQuadrant.c ----------- */

int getQuadrant(double mux, double muy, enum direction wind)
{
  /* --- Determine which quadrant ray points to --     -------------- */

  if (wind == DOWNWIND) {
    mux = -mux;
    muy = -muy;
  }

  if (mux >= 0.0) {
    if (muy >= 0.0) 
      return 1;  
    else 
      return 4; 
  } else {
    if (muy >= 0.0) 
      return 2; 
    else 
      return 3;
  }
}
/* ------- end ---------------------------- getQuadrant.c ----------- */

/* ------- begin -------------------------- getIntersect.c ---------- */

enum intersect getIntersect(double rx, double ry, double rz)
{
  enum intersect plane;
  double rmax;

  /* --- Determine which plane the ray intersects

   Note: r[xyz] = mu[xyz] / delta[xyz] --              -------------- */

  rx = fabs(rx);
  ry = fabs(ry);
  rz = fabs(rz);

  if (rx > ry) {
    plane = YZ;
    rmax  = rx;
  } else {
    plane = XZ;
    rmax  = ry;
  }
  if (rz > rmax) plane = XY;
 
  return plane;  
}
/* ------- end ------------------------------- getIntersect.c ------- */

/* ------- begin ----------------------------- calcFrac.c ----------- */

void calcFrac(double *frac1, double *frac2,
	      double rx, double ry, double rz, enum intersect plane,
	      enum direction wind)
{
  /* --- Calculate coordinate of point where the ray hits cell 
         boundary in units of cell size, relative to point ray 
         originates from.

    Note: -1 < frac[12] < 1
     --                                                -------------- */
 
  switch (plane) {
  case XY:
    *frac1 = rx/fabs(rz);
    *frac2 = ry/fabs(rz);
    break;
  case XZ:
    *frac1 = rx/fabs(ry);
    *frac2 = rz/fabs(ry);
    break;
  case YZ:
    *frac1 = ry/fabs(rx);
    *frac2 = rz/fabs(rx);
  }
  if (wind == DOWNWIND) {
    *frac1 *= -1.0;
    *frac2 *= -1.0;
  }
}
/* ------- end ---------------------------- calcFrac.c -------------- */

/* ------- begin -------------------------- calcDistance.c ---------- */

double calcDistance(double frac1, double frac2, enum intersect plane,
		    double dx, double dy, double dz)
{
  /* --- Calculate distance between point where ray originates from and 
         point where ray intersects nearest plane --   -------------- */

  switch (plane){
  case XY:
    return sqrt(SQ(dz) + SQ(frac1*dx) + SQ(frac2*dy));
  case XZ:
    return sqrt(SQ(dy) + SQ(frac1*dx) + SQ(frac2*dz));
  case YZ:
    return sqrt(SQ(dx) + SQ(frac1*dy) + SQ(frac2*dz));  
  }
}
/* ------- end ---------------------------- calcDistance.c ---------- */

/* ------- begin -------------------------- setStencil.c ------------ */

void setStencil(Geometry *geometry, int mu, int k, enum direction wind)
{
  const char routineName[] = "setStencil";
  register int l, m, n;

  int     l0, m0, k0, countXZ, countYZ, count, Nvalid;
  bool_t *valid;
  double  rx, ry, rz, frac1, frac2, muz, dz, longfracx, longfracy,
          svalid;
  Stencil *st;
  Longchar *lc;
  Intersect *vis;
  
  /* --- Fill interpolation stencil:

         Determine largest grid coordinate with x1_grid < x1_intersect 
         and x2_grid < x2_intersect, then fill stencil with addresses 
         of interpolation points surrounding the point of intersection.
     --                                                -------------- */

  if (wind == UPWIND) {
    st = &geometry->stencil_uw[mu][k];
    dz = geometry->z[k - 1] - geometry->z[k];
  } else {
    st = &geometry->stencil_dw[mu][k];
    dz = geometry->z[k] - geometry->z[k + 1];
  }
  muz = sqrt(1.0 - (SQ(geometry->mux[mu]) + SQ(geometry->muy[mu])));
  rx  = geometry->mux[mu] / geometry->dx;
  ry  = geometry->muy[mu] / geometry->dy;
  rz  = muz / dz;

  /* --- Calculate quadrant ray point to --            -------------- */
  
  st->quadrant = getQuadrant(geometry->mux[mu], geometry->muy[mu], wind);

  /* --- Determine which is the nearest plane in direction of ray - - */

  st->plane = getIntersect(rx, ry, rz);

  /* --- Determine fractions and distance between gridpoint and
         point where ray intersects nearest plane --   -------------- */

  calcFrac(&frac1, &frac2, rx, ry, rz, st->plane, wind);
  st->ds = calcDistance(frac1, frac2, st->plane,
			geometry->dx, geometry->dy, dz);

  switch (st->plane){
  case XY:
    st->zbase[0] = (wind == UPWIND) ?  k - 1  :  k + 1;

    l0 = floor(frac1);
    frac1 -= l0;

    m0 = floor(frac2);
    frac2 -= m0;

    switch (input.interpolate_3D) {
    case LINEAR_3D:
      st->xbase[0] = l0;
      st->xbase[1] = l0 + 1;
      st->xkernel[0] = 1.0 - frac1;
      st->xkernel[1] = frac1;

      st->ybase[0] = m0;
      st->ybase[1] = m0 + 1;
      st->ykernel[0] = 1.0 - frac2;
      st->ykernel[1] = frac2;
      break;

    case BICUBIC_3D:
      for (l = 0;  l < NCC;  l++) st->xbase[l] = l0 - 1 + l;
      cc_kernel(frac1, st->xkernel);

      for (m = 0;  m < NCC;  m++) st->ybase[m] = m0 - 1 + m;
      cc_kernel(frac2, st->ykernel);
    }
    break;

  case XZ:
    st->ybase[0] = (st->quadrant == 1 || st->quadrant == 2) ?  1  :  -1;

    l0 = floor(frac1);
    frac1 -= l0;

    if (frac2 < 0.0) {
      k0 = k + 1;
      frac2 += 1.0;
    } else
      k0 = k;

    st->xbase[0] = l0;
    st->xbase[1] = l0 + 1;
    st->xkernel[0] = 1.0 - frac1;
    st->xkernel[1] = frac1;

    st->zbase[0] = k0;
    st->zbase[1] = k0 - 1;
    st->zkernel[0] = 1.0 - frac2;
    st->zkernel[1] = frac2;
    break;

  case YZ:
    st->xbase[0] = (st->quadrant == 1 || st->quadrant == 4) ?  1  :  -1;

    m0 = floor(frac1);
    frac1 -= m0;

    if (frac2 < 0.0) {
      k0 = k + 1;
      frac2 += 1.0;
    } else
      k0 = k;

    st->ybase[0] = m0;
    st->ybase[1] = m0 + 1;
    st->ykernel[0] = 1.0 - frac1;
    st->ykernel[1] = frac1;

    st->zbase[0] = k0;
    st->zbase[1] = k0 - 1;
    st->zkernel[0] = 1.0 - frac2;
    st->zkernel[1] = frac2;
  }

  st->longchar = NULL;

  /* --- If XZ or YZ plane is hit first, fill long characteristic  
         structures --                                 -------------- */

  if (st->plane == YZ || st->plane == XZ) {
    if (wind == UPWIND) {
      longfracx = rx / rz;
      longfracy = ry / rz;
    } else {
      longfracx = -rx / rz;
      longfracy = -ry / rz;
    }
    /* --- Determine number of crossings of YZ and XZ planes,
           respectively, before XY plane is crossed -- -------------- */

    countYZ = (int) floor(fabs(longfracx));
    countXZ = (int) floor(fabs(longfracy));

    /* --- Tabulate the vertical intersects --         -------------- */

    vis = (Intersect *) malloc((countYZ + countXZ) * sizeof(Intersect));

    count = 0;
    for (l = 0;  l < countYZ;  l++) {
      vis[count].plane = YZ;
      vis[count].s     = ((l+1) * geometry->dx) / fabs(geometry->mux[mu]);
      vis[count].x     = geometry->mux[mu] * vis[count].s;
      vis[count].y     = geometry->muy[mu] * vis[count].s;
      vis[count].z     = muz * vis[count].s;
      count++;
    }
    for (m = 0;  m < countXZ;  m++) {
      vis[count].plane = XZ;
      vis[count].s     = ((m+1) * geometry->dy) / fabs(geometry->muy[mu]);
      vis[count].x     = geometry->mux[mu] * vis[count].s;
      vis[count].y     = geometry->muy[mu] * vis[count].s;
      vis[count].z     = muz * vis[count].s;
      count++;
    }
    /* --- Sort the intersections according to distance s from the
           point of origin --                          -------------- */

    qsort(vis, count, sizeof(Intersect), sascend);

    /* --- Determine if any of the intersections appear too close to
           one another, for instance when they occur on an intersection
           between YZ and XZ planes --                 -------------- */

    Nvalid = count;
    valid  = (bool_t *) malloc(Nvalid * sizeof(bool_t));
    valid[0] = TRUE;
    for (n = 1;  n < count;  n++) {
      if ((vis[n].s - vis[n-1].s) <
	  MIN_FRAC * MIN(geometry->dx, geometry->dy)) {
	valid[n] = FALSE;
	Nvalid--;
      } else
	valid[n] = TRUE;
    }

    st->longchar = (Longchar *) malloc(sizeof(Longchar));
    lc = st->longchar;
    lc->Nst = Nvalid + 1;
    lc->stencil = (Stencil *) malloc(lc->Nst * sizeof(Stencil));

    /* --- The first stencil to be stored is the intersection with
           the horizontal XY plane, ds is the distance to the first
           valid vertical intersection --              -------------- */

    lc->stencil[0].plane = XY;
    n = count;
    do {
      n--;
      lc->stencil[0].ds = dz / muz - vis[n].s;
    } while (!valid[n]);
    svalid = vis[n].s;

    l0 = floor(longfracx);
    frac1 = longfracx - l0;
    m0 = floor(longfracy);
    frac2 = longfracy - m0;

    switch (input.interpolate_3D) {
    case LINEAR_3D:
      lc->stencil[0].xbase[0] = l0;
      lc->stencil[0].xbase[1] = l0 + 1;
      lc->stencil[0].xkernel[0] = 1.0 - frac1;
      lc->stencil[0].xkernel[1] = frac1;

      lc->stencil[0].ybase[0] = m0;
      lc->stencil[0].ybase[1] = m0 + 1;
      lc->stencil[0].ykernel[0] = 1.0 - frac2;
      lc->stencil[0].ykernel[1] = frac2;

      break;
    case BICUBIC_3D:
      for (l = 0;  l < NCC;  l++) lc->stencil[0].xbase[l] = l0 - 1 + l;
      cc_kernel(frac1, lc->stencil[0].xkernel);

      for (m = 0;  m < NCC;  m++) lc->stencil[0].ybase[m] = m0 - 1 + m;
      cc_kernel(frac2, lc->stencil[0].ykernel);
    }
    lc->stencil[0].zbase[0] = (wind == UPWIND) ?  k - 1  :  k + 1;

    count = 1;
    for (n = (countXZ + countYZ)-1;  n >= 0;  n--) {
      if (valid[n]) {
	lc->stencil[count].plane = vis[n].plane;
	if (count == Nvalid)
	  lc->stencil[count].ds = svalid;
	else {
	  lc->stencil[count].ds = svalid - vis[n-1].s;
	  svalid = vis[n-1].s;
	}

	switch (vis[n].plane) {
        case XZ:
          l0 = floor(vis[n].x / geometry->dx);
          lc->stencil[count].xbase[0] = l0;
          lc->stencil[count].xbase[1] = l0 + 1;
          frac1 = vis[n].x / geometry->dx - l0;
          lc->stencil[count].xkernel[0] = 1.0 - frac1;
	  lc->stencil[count].xkernel[1] = frac1;

          lc->stencil[count].ybase[0] =
	    (int) round(vis[n].y / geometry->dy);
	  break;

        case YZ:
          lc->stencil[count].xbase[0] =
	    (int) round(vis[n].x / geometry->dx);

          m0 = floor(vis[n].y / geometry->dy);
          lc->stencil[count].ybase[0] = m0;
          lc->stencil[count].ybase[1] = m0 + 1;
          frac1 = vis[n].y / geometry->dy - m0;
          lc->stencil[count].ykernel[0] = 1.0 - frac1;
	  lc->stencil[count].ykernel[1] = frac1;
          break;

	case XY:
	  sprintf(messageStr, "Plane not implemented: %d", vis[n].plane);
	  Error(ERROR_LEVEL_2, routineName, messageStr);
	}

	if (wind == UPWIND) {
	  lc->stencil[count].zbase[0] = k;
	  lc->stencil[count].zbase[1] = k - 1;
	  frac2 = vis[n].z / dz;
	} else {
	  lc->stencil[count].zbase[0] = k + 1;
	  lc->stencil[count].zbase[1] = k;
	  frac2 = 1.0 - vis[n].z / dz;
	}
	lc->stencil[count].zkernel[0] = 1.0 - frac2;
	lc->stencil[count].zkernel[1] = frac2;

	count++;
      }
    }
    /* --- Free temporary arrays --                    -------------- */

    free(vis);
    free(valid);

    /* --- Store total number of long characteristics to be used - -- */

    if (st->plane == YZ)
      geometry->Nlongchar += geometry->Ny;
    else
      geometry->Nlongchar += geometry->Nx;

  }
}
/* ------- end ---------------------------- setStencil.c ------------ */

/* ------- begin -------------------------- round.c ----------------- */

double round(double x)
{
  double integr;
  
  modf(x + 0.5, &integr);
  return integr;
}
/* ------- end ---------------------------- round.c ----------------- */

/* ------- begin -------------------------- sascend.c --------------- */

int sascend(const void *vis1, const void *vis2)
{
  double s1 = ((Intersect *) vis1)->s, s2 = ((Intersect *) vis2)->s;
  
  if (s1 > s2)
    return 1;
  else if (s1 == s2)
    return 0;
  else
    return -1;
}
/* ------- end ---------------------------- sascend.c --------------- */

/* ------- begin -------------------------- Interpolate_3D.c -------- */

double Interpolate_3D(double *f, Geometry *geometry, Stencil *st,
		      int l, int m)
{
  register int i, j;

  int    index, offset1, offset2, xbase[NCC], ybase[NCC], zbase[2];
  bool_t monotonic = TRUE;
  double value;

  value = 0.0;

  switch (st->plane) {
  case XY:
    offset1 = st->zbase[0] * geometry->Nplane;

    switch (input.interpolate_3D) {
    case LINEAR_3D:
      for (i = 0;  i < 2;  i++) {
	xbase[i] = MODULO(l + st->xbase[i], geometry->Nx);
	ybase[i] = MODULO(m + st->ybase[i], geometry->Ny);
      }
      for (j = 0;  j < 2;  j++) {
	offset2 = offset1 + ybase[j] * geometry->Nx;
	for (i = 0;  i < 2;  i++)
	  value += f[offset2 + xbase[i]] *
	    st->xkernel[i] * st->ykernel[j];
      }
      break;
	   
    case BICUBIC_3D:
      for (i = 0;  i < NCC;  i++) {
	xbase[i] = MODULO(l + st->xbase[i], geometry->Nx);
	ybase[i] = MODULO(m + st->ybase[i], geometry->Ny);
      }
      value = cubicconvol(f + offset1, geometry->Nx, xbase, ybase,
			  st->xkernel, st->ykernel, monotonic);
    }
    break;

  case XZ:
    for (i = 0;  i < 2;  i++)
      xbase[i] = MODULO(l + st->xbase[i], geometry->Nx);
    offset1 = MODULO(m + st->ybase[0], geometry->Ny) * geometry->Nx;

    for (j = 0;  j < 2;  j++) {
      offset2 = offset1 + st->zbase[j] * geometry->Nplane;
      for (i = 0;  i < 2;  i++)
	value += f[offset2 + xbase[i]] * st->xkernel[i] * st->zkernel[j];
    }
    break;

  case YZ:
    offset1 = MODULO(l + st->xbase[0], geometry->Nx);
    for (i = 0;  i < 2;  i++)
      ybase[i] = MODULO(m + st->ybase[i], geometry->Ny);

    for (j = 0;  j < 2;  j++) {
      offset2 = st->zbase[j] * geometry->Nplane + offset1;
      for (i = 0;  i < 2;  i++) {
	index  = offset2 + ybase[i] * geometry->Nx;
	value += f[index] * st->ykernel[i] * st->zkernel[j];
      }
    }
  }  

  return value;
}
/* ------- end ---------------------------- Interpolate_3D.c -------- */

/* ------- begin -------------------------- cubicconvol.c ----------- */

double cubicconvol(double *f, int Nx, int *x_base, int *y_base,
                   double *x_kernel, double *y_kernel, bool_t monotonic)
{
  register int l, m, n;

  int    index;
  double result = 0.0, result_x, subA[NCC*NCC], A_min, A_max;

  /* --- Interpolation with cubic convolution on the subgrid
         x_base # y_base. If keyword monotonic is set the routine
         makes sure the interpolation value is between the min and max
         values of the inner part of the subgrid --    ------------- */

  for (m = 0, l = 0;  m < NCC;  m++) {
    result_x = 0.0;
    for (n = 0;  n < NCC;  n++, l++) {
      index = y_base[m]*Nx + x_base[n];
      subA[l] = f[index];
      result_x += x_kernel[n] * f[index];
    }
    result += y_kernel[m] * result_x;
  }
  if (monotonic) {
    A_min = MIN(MIN(subA[5], subA[6]), MIN(subA[9], subA[10]));
    A_max = MAX(MAX(subA[5], subA[6]), MAX(subA[9], subA[10]));

    if (result > A_max)
      return A_max;
    else if (result < A_min)
      return A_min;
    else
      return result;
  } else
    return result;
}
/* ------- end ---------------------------- cubicconvol.c ----------- */
