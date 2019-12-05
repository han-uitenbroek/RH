/* ------- file: -------------------------- piecewise_2D.c ----------

       Version:       rh2.0, 2-D Cartesian
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Wed Jun  6 14:10:11 2018 --

       --------------------------                      ----------RH-- */

/* --- Formal 2-D radiative transfer solver using the method of short
       characteristics. The transfer equation is solved in one direction
       along the ray only. Use keyword to_observer=TRUE to solve in the
       upward direction, and to_observer=FALSE to solve in opposite
       direction.

  See: P. B. Kunasz and  L. H. Auer 1988, JQSRT 39, 67-79

       For now it implements fixed and periodic boundary
       conditions in the horizontal direction.

  See: L. H. Auer, P. Fabiani Bendicho and J. Trujillo Bueno 1994,
       A&A, 292, 599-615

       At the top and bottom only a fixed boundary is allowed.

  See: L. H. Auer 2003, Formal Solution: EXPLICIT Answers, in 
       Stellar Atmosphere Modeling, ASP Conference Proceedings, Vol. 288,
       eds. Hubeny, Mihalas & Werner.


 Note: Rays with to_observer == TRUE go from BOTTOM to TOP.

       --                                              -------------- */

 
#include <math.h>
#include <stdlib.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "spectrum.h"
#include "geometry.h"
#include "inputs.h"
#include "error.h"
#include "bezier.h"


/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern Spectrum spectrum;
extern InputData input;


/* ------- begin -------------------------- Quadr.c ----------------- */

/* --- Evaluates function value with polynomial interpolation. -- --- */

double Quadr(double *values, Stencil *stencil,
	     enum direction direction, bool_t monotonic)
{
  const char routineName[] = "Quadr";
  register int n;

  int    *index;
  double  result, v[3], v_max, v_min, fraction, *coeff; 

  /* --- If monotonic = TRUE the interpolated value is forced to be
         monotonic with the original values.

    See: Auer & Paletou, A&A, 285, 675-686 --          -------------- */

  index    = stencil->index[direction];
  coeff    = stencil->coeff[direction];
  fraction = stencil->fraction[direction];

  result = 0.0;
  switch (stencil->order[direction]) {
  case SC_LINEAR:

    result = (1.0 - fraction)*values[index[0]] + fraction*values[index[1]];
    break;

  case SC_QUADRATIC:

    if (monotonic) {
      for (n = 0;  n < 3;  n++) {
	v[n]    = values[index[n]];
	result += coeff[n]*v[n];
      }
      v_max = MAX(MAX(v[0], v[1]), v[2]);
      v_min = MIN(MIN(v[0], v[1]), v[2]);

      if (result > v_max || result < v_min) {
	if (stencil->triplet[direction] == UPPER_TRIPLET)
	  result =
	    (1.0 - fraction)*values[index[0]] + fraction*values[index[1]];
	else
	  result =
	    (1.0 - fraction)*values[index[1]] + fraction*values[index[2]];
      }
    } else {
      for (n = 0;  n < 3;  n++)
	result += coeff[n] * values[index[n]];
    }
    break;

  default:
    Error(ERROR_LEVEL_2, routineName, "Illegal interpolation order");
  }

  return result;
}
/* ------- end ---------------------------- Quadr.c ----------------- */


/* ------- begin -------------------------- Piecewise_2D.c ---------- */

void Piecewise_2D(Geometry *geometry, int nspect, int mu,
		  bool_t to_observer, double *chi, double *S,
		  double *I, double *Psi)
{
  /* --- Parabolic interpolation of source function. -- ------------- */
  
  const char routineName[] = "Piecewise_2D";
  register int l, k, n;

  bool_t  monotonic = TRUE;
  enum    boundval bvalue;
  enum    direction upwind, downwind;
  enum    sweep sweep;
  int     Nx = geometry->Nx, Nz = geometry->Nz, local, dl, lstart,
          Nrow, row_left, row_right, dk, kstart, kend;
  double  mux = geometry->mux[mu], I_uw, *Iboundary, dtau_uw, dtau_dw,
          c1, c2, dS_uw, dS_dw, S_uw, chi_uw, S_dw, chi_dw, w[3],
          Sc, Bnu_uw, Bnu, T_uw, psi_uw, psi_dw, psi_0;
  Stencil  *st;
  LongChar *lc = NULL;

  /* --- In case of FIXED left and right boundary conditions treat
         the left and right columns --                 -------------- */

  if (geometry->hboundary == FIXED) {
    if ((to_observer && mux > 0.0) || (!to_observer && mux < 0.0)) {
      
      /* --- Rays that begin at the LEFT boundary and end at the RIGHT
	     boundary --                               -------------- */

      for (k = 0;  k < Nz;  k++) I[k*Nx] = geometry->Ileft[nspect][k];
      if (Psi) for (k = 0;  k < Nz;  k++) Psi[k*Nx] = 0.0;
    } else if ((to_observer && mux < 0.0) || (!to_observer && mux > 0.0)) {

      /* --- Rays that begin at the RIGHT boundary and end at the
             LEFT boundary --                          -------------- */ 

      for (k = 0;  k < Nz;  k++)
	I[(k+1)*Nx - 1] = geometry->Iright[nspect][k];
      if (Psi) for (k = 0;  k < Nz;  k++) Psi[(k+1)*Nx - 1] = 0.0;
    }
  }
  /* --- Boundary conditions in row at the top or bottom -- --------- */

  if (to_observer) {
    row_left = (Nz-1)*Nx;
    bvalue = geometry->bvalue[BOTTOM];
  } else {
    row_left = 0;
    bvalue = geometry->bvalue[TOP];
  }
  row_right = row_left + Nx;

  switch (bvalue) {
  case IRRADIATED:
    if (to_observer)
      Iboundary = geometry->Ibottom[nspect];
    else
      Iboundary = geometry->Itop[nspect];
    for (l = 0;  l < Nx;  l++)  I[row_left + l] = Iboundary[l];
    break;
  case ZERO:
    for (l = row_left;  l < row_right;  l++) I[l] = 0.0;
    break;
  case THERMALIZED:
    if (to_observer) {
      local = (Nz-1) * Nx;
      for (l = 0;  l < Nx;  l++) {
	st = &geometry->stencil[mu][local];

        T_uw = Quadr(atmos.T, st, UPWIND, monotonic);    
	Planck(1, &T_uw, spectrum.lambda[nspect], &Bnu_uw);
	Planck(1, &atmos.T[local], spectrum.lambda[nspect], &Bnu);

	chi_uw  = Quadr(chi, st, UPWIND, monotonic);
	dtau_uw = 0.5 * (chi_uw + chi[local]) * st->ds[UPWIND];

	I[local] = Bnu - (Bnu_uw - Bnu) / dtau_uw;
        local++;
      }
    } else {
      Error(ERROR_LEVEL_2, routineName, 
	    "Boundary condition THERMALIZED not implemented for TOP");
    }
    break;
  }

  if (Psi) for (l = row_left;  l < row_right;  l++) Psi[l] = 0.0;

  /* --- Finally, go through the grid to solve the non-local transfer
         equation --                                   -------------- */

  if (to_observer) {
    if (mux >= 0) dl = 1;  else  dl= -1;
    sweep  = UP;
    dk     = -1;
    kstart = Nz - 2;
    kend   = 0;
    upwind   = DOWNWIND;
    downwind = UPWIND;
  } else {
    if (mux <= 0) dl = 1;  else  dl= -1;
    sweep  = DOWN;
    dk     = 1;
    kstart = 1;
    kend   = Nz - 1;
    upwind   = UPWIND;
    downwind = DOWNWIND;
  }

  for (k = kstart;  k != kend + dk;  k += dk) {

    if (geometry->hboundary == PERIODIC) {
      lc = &geometry->longchar[mu][sweep*Nz + k];
      lstart = lc->lstart;
      Nrow = Nx;
    } else {
      if (mux == 0.0) {
        lstart = 0;
	Nrow = Nx;
      } else {
	if ((to_observer && mux > 0.0) || (!to_observer && mux < 0.0))
	  lstart = 1;
	else
	  lstart = Nx - 2;

        Nrow = Nx - 1;
      }
    }

    for (n = 0, l = lstart;  n < Nrow;  n++, l += dl) {
      if (geometry->hboundary == PERIODIC)
	local = k*Nx + MODULO(l, Nx);
      else
	local = k*Nx + l;
      st = &geometry->stencil[mu][local];

      chi_uw  = Quadr(chi, st, upwind, monotonic);
      S_uw    = Quadr(S, st, upwind, monotonic);
      dtau_uw = 0.5 * (chi_uw + chi[local]) * st->ds[upwind];
      dS_uw   = (S_uw - S[local]) / dtau_uw;

      if (geometry->hboundary == PERIODIC &&
	  l == lc->lstart && lc->Nlc > 0)
	I_uw = SolveLong(lc, local, chi, S, I);
      else
	I_uw = Quadr(I, st, upwind, monotonic);
      
      if (k == kend ||
	  (geometry->hboundary == FIXED && n == Nrow-1)) {
	w2(dtau_uw, w);

	/* --- Piecewise linear at end of ray --     -------------- */

	I[local] = (1.0 - w[0])*I_uw + w[0]*S[local] + w[1]*dS_uw;
	if (Psi) Psi[local] = w[0] - w[1]/dtau_uw;

      } else {
	w3(dtau_uw, w);

	/* --- Piecewise quadratic elsewhere --      -------------- */

	chi_dw  = Quadr(chi, st, downwind, monotonic);
	dtau_dw = 0.5 * (chi[local] + chi_dw) * st->ds[downwind];
	S_dw    = Quadr(S, st, downwind, monotonic);
	dS_dw   = (S[local] - S_dw) / dtau_dw;

	c1 = (dS_uw*dtau_dw + dS_dw*dtau_uw);
	c2 = dS_uw - dS_dw;

	I[local] = (1.0 - w[0])*I_uw + w[0]*S[local] +
	  (w[1]*c1 + w[2]*c2) / (dtau_uw + dtau_dw);

	/* --- Try piecewise linear if quadratic gives negative
	       monochromatic intensity --              ------------ */ 

	if (I[local] < 0.0) {
	  I[local] = (1.0 - w[0])*I_uw + w[0]*S[local] + w[1]*dS_uw;

	  if (Psi) Psi[local] = w[0] - w[1]/dtau_uw;
	} else {
	  if (Psi) {
	    c1 = dtau_uw - dtau_dw;
	    Psi[local] = w[0] + (w[1]*c1 - w[2]) / (dtau_uw * dtau_dw);
	  }
	}
      }
    }
  }
}
/* ------- end ---------------------------- Piecewise_2D.c ---------- */


/* ------- begin -------------------------- Piecewise_Linear_2D.c --- */

void Piecewise_Linear_2D(Geometry *geometry, int nspect, int mu,
			 bool_t to_observer, double *chi, double *S,
			 double *I, double *Psi)
{
  /* --- Linear interpolation of source function --    -------------- */
 
  const char routineName[] = "Piecewise_Linear_2D";
  register int l, k, n;

  bool_t  monotonic = TRUE;
  enum    boundval bvalue;
  enum    direction upwind, downwind;
  enum    sweep sweep;
  int     Nx = geometry->Nx, Nz = geometry->Nz, local, dl, lstart,
          Nrow, row_left, row_right, dk, kstart, kend;
  double  mux = geometry->mux[mu], I_uw, *Iboundary, dtau_uw, dtau_dw,
          c1, dS_dw, S_uw, chi_uw, S_dw, chi_dw, w[2],
          Sc, Bnu_uw, Bnu, T_uw, psi_uw, psi_dw, psi_0;
  Stencil  *st;
  LongChar *lc = NULL;

  /* --- In case of FIXED left and right boundary conditions treat
         the left and right columns --                 -------------- */

  if (geometry->hboundary == FIXED) {
    if ((to_observer && mux > 0.0) || (!to_observer && mux < 0.0)) {
      
      /* --- Rays that begin at the LEFT boundary and end at the RIGHT
	     boundary --                               -------------- */

      for (k = 0;  k < Nz;  k++) I[k*Nx] = geometry->Ileft[nspect][k];
      if (Psi) for (k = 0;  k < Nz;  k++) Psi[k*Nx] = 0.0;
    } else if ((to_observer && mux < 0.0) || (!to_observer && mux > 0.0)) {

      /* --- Rays that begin at the RIGHT boundary and end at the
             LEFT boundary --                          -------------- */ 

      for (k = 0;  k < Nz;  k++)
	I[(k+1)*Nx - 1] = geometry->Iright[nspect][k];
      if (Psi) for (k = 0;  k < Nz;  k++) Psi[(k+1)*Nx - 1] = 0.0;
    }
  }
  /* --- Boundary conditions in row at the top or bottom -- --------- */

  if (to_observer) {
    row_left = (Nz-1)*Nx;
    bvalue = geometry->bvalue[BOTTOM];
  } else {
    row_left = 0;
    bvalue = geometry->bvalue[TOP];
  }
  row_right = row_left + Nx;

  switch (bvalue) {
  case IRRADIATED:
    if (to_observer)
      Iboundary = geometry->Ibottom[nspect];
    else
      Iboundary = geometry->Itop[nspect];
    for (l = 0;  l < Nx;  l++)  I[row_left + l] = Iboundary[l];
    break;
  case ZERO:
    for (l = row_left;  l < row_right;  l++) I[l] = 0.0;
    break;
  case THERMALIZED:
    if (to_observer) {
      local = (Nz-1) * Nx;
      for (l = 0;  l < Nx;  l++) {
	st = &geometry->stencil[mu][local];

        T_uw = Quadr(atmos.T, st, UPWIND, monotonic);    
	Planck(1, &T_uw, spectrum.lambda[nspect], &Bnu_uw);
	Planck(1, &atmos.T[local], spectrum.lambda[nspect], &Bnu);

	chi_uw  = Quadr(chi, st, UPWIND, monotonic);
	dtau_uw = 0.5 * (chi_uw + chi[local]) * st->ds[UPWIND];

	I[local] = Bnu - (Bnu_uw - Bnu) / dtau_uw;
        local++;
      }
    } else {
      Error(ERROR_LEVEL_2, routineName, 
	    "Boundary condition THERMALIZED not implemented for TOP");
    }
    break;
  }

  if (Psi) for (l = row_left;  l < row_right;  l++) Psi[l] = 0.0;

  /* --- Finally, go through the grid to solve the non-local transfer
         equation --                                   -------------- */

  if (to_observer) {
    if (mux >= 0) dl = 1;  else  dl= -1;
    sweep  = UP;
    dk     = -1;
    kstart = Nz - 2;
    kend   = 0;
    upwind   = DOWNWIND;
    downwind = UPWIND;
  } else {
    if (mux <= 0) dl = 1;  else  dl= -1;
    sweep  = DOWN;
    dk     = 1;
    kstart = 1;
    kend   = Nz - 1;
    upwind   = UPWIND;
    downwind = DOWNWIND;
  }

  for (k = kstart;  k != kend + dk;  k += dk) {

    if (geometry->hboundary == PERIODIC) {
      lc = &geometry->longchar[mu][sweep*Nz + k];
      lstart = lc->lstart;
      Nrow = Nx;
    } else {
      if (mux == 0.0) {
        lstart = 0;
	Nrow = Nx;
      } else {
	if ((to_observer && mux > 0.0) || (!to_observer && mux < 0.0))
	  lstart = 1;
	else
	  lstart = Nx - 2;

        Nrow = Nx - 1;
      }
    }

    for (n = 0, l = lstart;  n < Nrow;  n++, l += dl) {
      if (geometry->hboundary == PERIODIC)
	local = k*Nx + MODULO(l, Nx);
      else
	local = k*Nx + l;
      st = &geometry->stencil[mu][local];

      chi_uw  = Quadr(chi, st, upwind, monotonic);
      dtau_uw = 0.5 * (chi_uw + chi[local]) * st->ds[upwind];     
      S_uw    = Quadr(S, st, upwind, monotonic);

      if (geometry->hboundary == PERIODIC && l == lc->lstart && lc->Nlc > 0)
	I_uw = SolveLong(lc, local, chi, S, I);
      else
	I_uw = Quadr(I, st, upwind, monotonic);
      
      w2(dtau_uw, w);

      c1 = (S_uw - S[local]) / dtau_uw;
      I[local] = (1.0 - w[0])*I_uw + w[0]*S[local] + w[1]*c1;

      if (Psi) Psi[local] = w[0] - w[1]/dtau_uw;
    }
  }
}
/* ------- end ---------------------------- Piecewise_Linear_2D.c --- */
