/* ------- file: -------------------------- piecestokes_3D.c --------

       Version:       rh2.0, 3-D Cartesian, short characteristics
       Author:        Han Uitenbroek(huitenbroek@nso.edu)
       Last modified: Tue Sep 23 12:29:54 2025 --

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


/* ------- begin -------------------------- ShortChar_Stokes.c ------ */

void ShortChar_Stokes(Geometry *geometry, int nspect, int mu,
		      bool_t to_observer,
		      double *chi, double **S, double **I, double *Psi)
{
  const char routineName[] = "ShortChar_Stokes";
  register int k, l, m, n;

  int     Nz = geometry->Nz, Nplane = geometry->Nplane, quadrant,
          dl, lstart, lend, dm, mstart, mend, dk, kstart, kend, local;
  double  I_uw[4], T_uw, Bnu_uw, Bnu, chi_uw, dtau_uw;
  Stencil *st_uw, *st_dw;

  if (to_observer) {

    /* --- Boundary condition at the BOTTOM of the grid -- ---------- */

    switch (geometry->z_boundary_value[BOTTOM]) {
    case IRRADIATED:
      for (l = 0;  l < Nplane;  l++)
	I[0][(Nz-1)*Nplane + l] = geometry->Ibottom[nspect][l];
      break;

    case ZERO:
      for (l = (Nz-1)*Nplane;  l < Nz*Nplane;  l++)
	I[0][l] = 0.0;
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

	    I[0][local] = Bnu - (Bnu_uw - Bnu) / dtau_uw;
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
    /* --- Assume irradiation is unpolarized --        -------------- */

    for (n = 1;  n < 4;  n++)
      for (l = (Nz-1)*Nplane;  l < Nz*Nplane;  l++) I[n][l] = 0.0;

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
      for (l = 0;  l < Nplane;  l++) I[0][l] = geometry->Itop[nspect][l];
      break;

    case ZERO:
      for (l = 0;  l < Nplane;  l++) I[0][l] = 0.0;
      break;

    case THERMALIZED:
      Planck(Nplane, atmos.T, spectrum.lambda[nspect], I[0]);
      break;

    case REFLECTIVE:
      Error(ERROR_LEVEL_2, routineName, 
	    "Boundary condition REFLECTIVE not implemented");
      break;
    }
    /* --- Assume irradiation is unpolarized --        -------------- */

    for (n = 1;  n < 4;  n++)
      for (l = 0;  l < Nplane;  l++) I[n][l] = 0.0;

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
	    SolveLongStokes(geometry, st_uw->longchar, nspect,
			    k, l, m, chi, S, I, I_uw);
	  else {
	    for (n = 0;  n < 4;  n++)
	      I_uw[n] = Interpolate_3D(I[n], geometry, st_uw, l, m);
	  }
	  if (input.S_interpolation_stokes == DELO_PARABOLIC) {
	    Piece_Stokes_3D(geometry, st_uw, st_dw, nspect,
			    k, kend, l, m, I_uw, chi, S, I, Psi);
	  } else if (input.S_interpolation_stokes == DELO_BEZIER3) {
	    Piece_Stokes_Bezier3_3D(geometry, st_uw, st_dw, nspect,
				    k, kend, l, m, I_uw, chi, S, I, Psi);
	  } else {
	    sprintf(messageStr,
		    "Unknown polarization solver: %d",
		    input.S_interpolation_stokes);
	    Error(ERROR_LEVEL_1, routineName, messageStr);
	  }
	}
      }
    } else {

      /* --- Inner loop is in the y-direction --       -------------- */      

      for (l = lstart;  l != lend+dl;  l += dl) {
	for (m = mstart;  m != mend+dm;  m += dm) {
	  if (l == lstart)
	    SolveLongStokes(geometry, st_uw->longchar, nspect,
			    k, l, m, chi, S, I, I_uw);
	  else {
	    for (n = 0;  n < 4;  n++)
	      I_uw[n] = Interpolate_3D(I[n], geometry, st_uw, l, m);
	  }
	  if (input.S_interpolation_stokes == DELO_PARABOLIC) {
	    Piece_Stokes_3D(geometry, st_uw, st_dw, nspect,
			    k, kend, l, m, I_uw, chi, S, I, Psi);
	  } else if (input.S_interpolation_stokes == DELO_BEZIER3) {
	    Piece_Stokes_Bezier3_3D(geometry, st_uw, st_dw, nspect,
				    k, kend, l, m, I_uw, chi, S, I, Psi);
	  } else {
	    sprintf(messageStr,
		    "Unknown polarization solver: %d",
		    input.S_interpolation_stokes);
	    Error(ERROR_LEVEL_1, routineName, messageStr);
	  }
	}
      }
    }
  }
}
/* ------- end ---------------------------- ShortChar_Stokes.c ------ */

/* ------- begin -------------------------- Piece_Stokes_3D.c ------- */

void Piece_Stokes_3D(Geometry *geometry, Stencil *st_uw, Stencil *st_dw,
		     int nspect, int k, int kend, int l, int m,
		     double *I_uw,
		     double *chi, double **S, double **I, double *Psi)
{
  /* --- Piecewise integration of the coupled Stokes transfer equations
         in two dimensions. Method is quasi-parabolic DELO method.

    See: - D. E. Rees, G. A. Murphy and C. J. Durrant 1989, ApJ 339,
           1093-1106.

         - H. Socas Navarro, J. Trujillo Bueno and B. Ruiz Cobo 2000,
           "Non-LTE Inversion of Stokes Profiles", ApJ 530, 977.
      --                                               -------------- */

  register int n, j;

  int    local;
  double chi_uw, chi_dw, S_uw[4], S_dw[4], dS_dw[4], dS_uw[4],
         dtau_uw, dtau_dw, w[3], c1, c2, P[4], Q[4][4], **R,
         K[4][4], K_uw[4][4]; 

  R = matrix_double(4, 4);

  local = k*geometry->Nplane + m*geometry->Nx + l;

  /* --- The upwind quantities --                      -------------- */

  chi_uw  = Interpolate_3D(chi, geometry, st_uw, l, m);
  dtau_uw = 0.5 * (chi_uw + chi[local]) * st_uw->ds;
  StokesK_3D(nspect, geometry, st_uw, l, m, chi_uw, K_uw);

  for (n = 0;  n < 4;  n++) {
    S_uw[n]  = Interpolate_3D(S[n], geometry, st_uw, l, m);
    dS_uw[n] = (S_uw[n] - S[n][local]) / dtau_uw;
  }
  StokesK(nspect, local, chi[local], K);

  if (k == kend) {
    w2(dtau_uw, w);

    /* --- Piecewise linear integration in last layer -- ------------ */

    for (n = 0;  n < 4;  n++) {
      c1 = (S_uw[n] - S[n][local]) / dtau_uw;
      P[n] = w[0]*S[n][local] + w[1]*dS_uw[n];
    }
    if (Psi) Psi[local] = w[0] - w[1]/dtau_uw;
  } else {
    w3(dtau_uw, w);

    /* --- The downwind quantities --                  -------------- */

    chi_dw  = Interpolate_3D(chi, geometry, st_dw, l, m);
    dtau_dw = 0.5 * (chi[local] + chi_dw) * st_dw->ds;

    /* --- Piecewise quadratic integration --          -------------- */

    for (n = 0;  n < 4;  n++) {
      S_dw[n]  = Interpolate_3D(S[n], geometry, st_dw, l, m);
      dS_dw[n] = (S[n][local] - S_dw[n]) / dtau_dw;
      c1 = dS_uw[n]*dtau_dw + dS_dw[n]*dtau_uw;
      c2 = dS_uw[n] - dS_dw[n];
      P[n] = w[0]*S[n][local] + (w[1]*c1 + w[2]*c2) /
	(dtau_uw + dtau_dw);
    }
    if (Psi) {
      c1 = dtau_uw - dtau_dw;
      Psi[local] = w[0] + (w[1]*c1 - w[2]) / (dtau_uw * dtau_dw);
    }
  }
  for (n = 0;  n < 4;  n++) {
    for (j = 0;  j < 4;  j++) {
      Q[n][j] = -w[1]/dtau_uw * K_uw[n][j];
      R[n][j] = (w[0] - w[1]/dtau_uw) * K[n][j];
    }
    Q[n][n] = 1.0 - w[0];
    R[n][n] = 1.0;
  }
  for (n = 0;  n < 4;  n++) {
    for (j = 0;  j < 4;  j++) 
      P[n] += Q[n][j] * I_uw[j];
  }
  /* --- Solve linear equations for I --               -------------- */
      
  SolveLinearEq(4, R, P, TRUE);
      
  /* --- Store results for Stokes vector --            -------------- */
      
  for (n = 0;  n < 4;  n++) I[n][local] = P[n];

  freeMatrix((void **) R);
}
/* ------- end ---------------------------- Piece_Stokes_3D.c ------- */

/* ------- begin -------------------------- StokesK_3D.c ------------ */

void StokesK_3D(int nspect, Geometry *geometry, Stencil *st,
		int l, int m, double chi_I, double K[4][4])
{
  register int i, j;

  ActiveSet *as;

  /* --- Return the elements of the reduced 4x4 Stokes
         opacity matrix K', which is defined as:

           =           =         =              =
           K' = (chi_c*1 + chi_l*Phi) / chi_I - 1,

	   for wavelength# nspect, spatial point k, and ray mu.

    See: Rees, Murphy, & Durrant, 1989, ApJ 339, 1093-1106.

         This is the 3-D short characteristics version.
         --                                            -------------- */
  as = &spectrum.as[nspect];

  for (j = 0;  j < 4;  j++)
    for (i = 0;  i < 4;  i++) K[j][i] = 0.0;

  /* --- First the contribution from the active set -- -------------- */

  if (containsPolarized(as)) {
    K[0][1] = Interpolate_3D(as->chi + atmos.Nspace,
                             geometry, st, l, m);
    K[0][2] = Interpolate_3D(as->chi + 2*atmos.Nspace,
                             geometry, st, l, m);
    K[0][3] = Interpolate_3D(as->chi + 3*atmos.Nspace,
                             geometry, st, l, m);

    if (input.magneto_optical) {
      K[1][2] = Interpolate_3D(as->chip + 2*atmos.Nspace,
                               geometry, st, l, m);
      K[1][3] = Interpolate_3D(as->chip + atmos.Nspace,
                               geometry, st, l, m);
      K[2][3] = Interpolate_3D(as->chip,
                               geometry, st, l, m);
    }
  }
  /* --- Add possible contrubution from the background -- ----------- */

  if (atmos.backgrflags[nspect].ispolarized) {
    K[0][1] += Interpolate_3D(as->chi_c + atmos.Nspace,
                              geometry, st, l, m);
    K[0][2] += Interpolate_3D(as->chi_c + 2*atmos.Nspace,
                              geometry, st, l, m);
    K[0][3] += Interpolate_3D(as->chi_c + 3*atmos.Nspace,
                              geometry, st, l, m);

    if (input.magneto_optical) {
      K[1][2] += Interpolate_3D(as->chip_c + 2*atmos.Nspace,
                                geometry, st, l, m);
      K[1][3] += Interpolate_3D(as->chip_c + atmos.Nspace,
                                geometry, st, l, m);
      K[2][3] += Interpolate_3D(as->chip_c,
                                geometry, st, l, m);
    }
  }
  /* --- Divide by Stokes I opacity and fill lower diagonal part -- - */

  for (j = 0;  j < 3;  j++) {
    for (i = j+1;  i < 4;  i++) {
      K[j][i] /= chi_I;
      K[i][j]  = K[j][i];
    }
  }
  /* --- Anti-symmetric magneto-optical elements --    -------------- */

  if (input.magneto_optical) {
    K[1][3] *= -1.0;
    K[2][1] *= -1.0;
    K[3][2] *= -1.0;
  }
}
/* ------- end ---------------------------- StokesK_3D.c ------------ */
