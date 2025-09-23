/* ------- file: -------------------------- piecewise_3D.c ----------

       Version:       rh2.0, 3-D Cartesian, short characteristics
       Author:        Han Uitenbroek(huitenbroek@nso.edu)
       Last modified: Tue Sep 23 12:28:26 2025 --

       --------------------------                      ----------RH-- */

#include <stdlib.h>
#include <math.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "geometry.h"
#include "spectrum.h"
#include "inputs.h"
#include "error.h"


/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern Spectrum spectrum;
extern InputData input;
extern char messageStr[];


/* ------- begin -------------------------- ShortChar.c ------------- */

void ShortChar(Geometry *geometry, int nspect, int mu,
               bool_t to_observer, double *chi, double *S,
               double *I, double *Psi)
{
  const char routineName[] = "ShortChar";
  register int k, l, m, ml;

  int     Nz = geometry->Nz, Nplane = geometry->Nplane, quadrant,
          dl, lstart, lend, dm, mstart, mend, dk, kstart, kend, local;
  double  I_uw, Bnu_uw, Bnu, T_uw, chi_uw, dtau_uw;
  Stencil *st_uw, *st_dw;

  if (to_observer) {

    /* --- Boundary condition at the BOTTOM of the grid -- ---------- */

    switch (geometry->z_boundary_value[BOTTOM]) {
    case IRRADIATED:
      for (l = 0;  l < Nplane;  l++)
	I[(Nz-1)*Nplane + l] = geometry->Ibottom[nspect][l];
      break;

    case ZERO:
      for (l = (Nz-1)*Nplane;  l < Nz*Nplane;  l++)
	I[l] = 0.0;
      break;

    case THERMALIZED:
      if (to_observer) {
        local = (geometry->Nz-1) * geometry->Nplane;
	st_uw = &geometry->stencil_uw[mu][geometry->Nz-1];

	for (m = 0;  m < geometry->Ny;  m++) {
	  for (l = 0;  l < geometry->Nx;  l++) {	    
	    T_uw = Interpolate_3D(atmos.T, geometry, st_uw, l, m);
	    Planck(1, &T_uw, spectrum.lambda[nspect], &Bnu_uw);
	    Planck(1, &atmos.T[local], spectrum.lambda[nspect], &Bnu);

	    chi_uw  = Interpolate_3D(chi, geometry, st_uw, l, m);
	    dtau_uw = 0.5 * (chi_uw + chi[local]) * st_uw->ds;

	    I[local] = Bnu - (Bnu_uw - Bnu) / dtau_uw;
	    local++;
	  }
	}
      } else {
	Error(ERROR_LEVEL_2, routineName, 
	      "Boundary condition THERMALIZED not implemented for TOP");
      }
      break;
    case REFLECTIVE:
      Error(ERROR_LEVEL_2, routineName, 
	    "Boundary condition REFLECTIVE not implemented");
      break;
    }
    if (Psi)
      for (l = (Nz-1)*Nplane;  l < Nz*Nplane;  l++) Psi[l] = 0.0;

    dk     = -1;
    kstart = Nz - 2;
    kend   = 0;

    quadrant = geometry->stencil_uw[mu][kstart].quadrant;
  } else {

    /* --- Boundary condition at the TOP of the grid -- ------------- */

    switch (geometry->z_boundary_value[TOP]) {
    case IRRADIATED:
      for (l = 0;  l < Nplane;  l++) I[l] = geometry->Itop[nspect][l];
      break;

    case ZERO:
      for (l = 0;  l < Nplane;  l++) I[l] = 0.0;
      break;

    case THERMALIZED:
      Planck(Nplane, atmos.T, spectrum.lambda[nspect], I);
      break;
      
    case REFLECTIVE:
      Error(ERROR_LEVEL_2, routineName, 
	    "Boundary condition REFLECTIVE not implemented");
      break;
    }
    if (Psi) for (l = 0;  l < Nplane;  l++) Psi[l] = 0.0;

    dk     = 1;
    kstart = 1;
    kend   = Nz - 1;

    quadrant = geometry->stencil_dw[mu][kstart].quadrant;
  }
  dl = (quadrant == 1 || quadrant == 4) ?  1  :  -1;
  dm = (quadrant == 1 || quadrant == 2) ?  1  :  -1;

  /* --- Figure out in which direction the grids have to be crossed - */

  if (dl == 1) {
    lstart = 0;
    lend   = geometry->Nx - 1;
  } else {
    lstart = geometry->Nx - 1;
    lend   = 0;
  }
  if (dm == 1) {
    mstart = 0;
    mend   = geometry->Ny - 1;
  } else {
    mstart = geometry->Ny - 1;
    mend   = 0;
  }

  for (k = kstart;  k != kend+dk;  k += dk) {
    if (to_observer) {
      st_uw = &geometry->stencil_dw[mu][k];
      st_dw = &geometry->stencil_uw[mu][k];
    } else {
      st_uw = &geometry->stencil_uw[mu][k];
      st_dw = &geometry->stencil_dw[mu][k];
    }
    if (st_uw->plane == XY || st_uw->plane == XZ) {

      /* --- Inner loop is in the x-direction --       -------------- */

      for (m = mstart;  m != mend+dm;  m += dm) {
	for (l = lstart;  l != lend+dl;  l += dl) {
	  if (st_uw->plane == XZ && m == mstart)
	    I_uw = SolveLong(geometry, st_uw->longchar,
			     k, l, m, chi, S, I);
	  else
	    I_uw = Interpolate_3D(I, geometry, st_uw, l, m);

	  if (input.S_interpolation == S_LINEAR) {
	    Piecewise_Linear_3D(geometry, st_uw, st_dw,
				k, kend, l, m, I_uw, chi, S, I, Psi);
	  } else if (input.S_interpolation == S_PARABOLIC) {
	    Piecewise_3D(geometry, st_uw, st_dw,
			 k, kend, l, m, I_uw, chi, S, I, Psi);
	  } else if (input.S_interpolation == S_BEZIER3) {
	    Piecewise_Bezier3_3D(geometry, st_uw, st_dw,
				k, kend, l, m, I_uw, chi, S, I, Psi);
	  } else {
	    sprintf(messageStr,
		    "Unknown radiation solver: %d",
		    input.S_interpolation);
	    Error(ERROR_LEVEL_1, routineName, messageStr);
	  }
	}
      }
    } else {

      /* --- Inner loop is in the y-direction --       -------------- */

      for (l = lstart;  l != lend+dl;  l += dl) {
	for (m = mstart;  m != mend+dm;  m += dm) {
	  if (l == lstart)
	    I_uw = SolveLong(geometry, st_uw->longchar,
			     k, l, m, chi, S, I);
	  else
	    I_uw = Interpolate_3D(I, geometry, st_uw, l, m);
	  
	  if (input.S_interpolation == S_LINEAR) {
	    Piecewise_Linear_3D(geometry, st_uw, st_dw,
				k, kend, l, m, I_uw, chi, S, I, Psi);
	  } else if (input.S_interpolation == S_PARABOLIC) {
	    Piecewise_3D(geometry, st_uw, st_dw,
			 k, kend, l, m, I_uw, chi, S, I, Psi);
	  } else if (input.S_interpolation == S_BEZIER3) {
	    Piecewise_Bezier3_3D(geometry, st_uw, st_dw,
				k, kend, l, m, I_uw, chi, S, I, Psi);
	  } else {
	    sprintf(messageStr,
		    "Unknown radiation solver: %d",
		    input.S_interpolation);
	    Error(ERROR_LEVEL_1, routineName, messageStr);
	  }
	}
      }
    }
  }
}
/* ------- end ---------------------------- ShortChar.c ------------- */

/* ------- begin -------------------------- Piecewise_3D.c ---------- */

void Piecewise_3D(Geometry *geometry, Stencil *st_uw, Stencil *st_dw,
		  int k, int kend, int l, int m, double I_uw,
		  double *chi, double *S, double *I, double *Psi)
{
  int    local;
  double chi_uw, chi_dw, S_uw, S_dw, dS_dw, dS_uw,
         dtau_uw, dtau_dw, w[3], c1, c2; 

  local = k*geometry->Nplane + m*geometry->Nx + l;

  /* --- The upwind quantities --                      -------------- */

  chi_uw  = Interpolate_3D(chi, geometry, st_uw, l, m);
  S_uw    = Interpolate_3D(S, geometry, st_uw, l, m);
  dtau_uw = 0.5 * (chi_uw + chi[local]) * st_uw->ds;
  dS_uw   = (S_uw - S[local]) / dtau_uw;

  if (k == kend) {
    w2(dtau_uw, w);

    /* --- Piecewise linear integration in last layer -- ------------ */

    I[local] = I_uw*(1.0 - w[0]) + w[0]*S[local] + w[1]*dS_uw;

    if (Psi) Psi[local] = w[0] - w[1]/dtau_uw;
  } else {
    w3(dtau_uw, w);

    /* --- The downwind quantities --                  -------------- */

    chi_dw  = Interpolate_3D(chi, geometry, st_dw, l, m);
    dtau_dw = 0.5 * (chi[local] + chi_dw) * st_dw->ds;
    S_dw    = Interpolate_3D(S, geometry, st_dw, l, m);
    dS_dw   = (S[local] - S_dw) / dtau_dw;

    /* --- Piecewise quadratic integration --          -------------- */

    c1 = (dS_uw*dtau_dw + dS_dw*dtau_uw);
    c2 = (dS_uw - dS_dw);

    I[local] = I_uw*(1.0 - w[0]) + w[0]*S[local] + 
      (w[1]*c1 + w[2]*c2) / (dtau_uw + dtau_dw);

    /* --- Try piecewise linear if quadratic gives negative
           monochromatic intensity --                  -------------- */ 

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
/* ------- end ---------------------------- Piecewise_3D.c ---------- */


/* ------- begin -------------------------- Piecewise_Linear_3D.c --- */

void Piecewise_Linear_3D(Geometry *geometry, Stencil *st_uw, Stencil *st_dw,
			 int k, int kend, int l, int m, double I_uw,
			 double *chi, double *S, double *I, double *Psi)
{
  int    local;
  double chi_uw, S_uw, dS_uw, dtau_uw, w[2], c1, c2; 

  local = k*geometry->Nplane + m*geometry->Nx + l;

  /* --- The upwind quantities --                      -------------- */

  chi_uw  = Interpolate_3D(chi, geometry, st_uw, l, m);
  S_uw    = Interpolate_3D(S, geometry, st_uw, l, m);
  dtau_uw = 0.5 * (chi_uw + chi[local]) * st_uw->ds;
  dS_uw   = (S_uw - S[local]) / dtau_uw;

  w2(dtau_uw, w);

  I[local] = (1.0 - w[0])*I_uw + w[0]*S[local] + w[1]*dS_uw;

  if (Psi) Psi[local] = w[0] - w[1]/dtau_uw;
}
/* ------- end ---------------------------- Piecewise_Linear_3D.c --- */
