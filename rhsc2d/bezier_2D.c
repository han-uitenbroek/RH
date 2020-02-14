/* ------- file: -------------------------- bezier_2D.c -------------

   Cubic DELO-Bezier (polarized) and cubic short-char Bezier solvers.
   
   References: de la Cruz Rodriguez & Piskunov (2013), Auer (2003)
               (Derivatives) Fritsch & Butland (1984),
	       
   1D version Coded by J. de la Cruz Rodriguez (ISP-SU 2017),
   adapted to 2D by Han Uitenbroek

       Last modified: Mon Jan  6 14:50:08 2020 --

       --------------------------                      ----------RH-- */


#include <math.h>
#include <string.h>

#include "rh.h"
#include "error.h"
#include "atom.h"
#include "atmos.h"
#include "geometry.h"
#include "spectrum.h"
#include "bezier.h"


/* --- Identity matrix --                              -------------- */

static const double ident[4][4] =
  {{1.0, 0.0, 0.0, 0.0},
   {0.0, 1.0, 0.0, 0.0},
   {0.0, 0.0, 1.0, 0.0},
   {0.0, 0.0, 0.0, 1.0}};


/* --- Global variables --                             -------------- */

extern Geometry geometry;
extern Atmosphere atmos;
extern Spectrum spectrum;
extern char messageStr[];


/* ------- begin -------------------------- Piece_Stokes_Bezier3_2D.c */

void Piece_Stokes_Bezier3_2D(Geometry *geometry, int nspect, int mu,
			     bool_t to_observer, double *chi, double **S,
			     double **I, double *Psi)
{
  /* --- Cubic DELO-Bezier solver for polarized light
         Coded by J. de la Cruz Rodriguez (ISP-SU 2017).
         Adapted to 2-D version by Han Uitenbroek.

         Reference(s):
         J. de la Cruz Rodriguez & N. Piskunov (2013)
         --                                        ------------------ */

  const char routineName[] = "Piecewise_Bezier3_2D";
  register int l, lp, k, n, m, i, j;
  
  bool_t  monotonic = TRUE;
  enum    boundval bvalue;
  enum    direction upwind, downwind;
  enum    sweep sweep;
  int     Nx = geometry->Nx, Nz = geometry->Nz;
  int     local, dl, lstart, row_left, row_right, dk, kstart, kend;
  double  mux = geometry->mux[mu], I_uw[4], dtau_uw, dtau_dw, c1, c2;
  double  chi_uw, chi_dw, w[2], Bnu_uw, Bnu, T_uw, *Iboundary;
  double  dchi_uw, dchi_c, dchi_dw, ds03, dt03;
  double  eps = 0, alpha = 0, beta = 0, gamma = 0, theta = 0;
  double  K_uw[4][4], K_loc[4][4], K_dw[4][4];
  double  dK_uw[4][4], dK_c[4][4];
  double  S_uw[4], S_loc[4], S_dw[4], dS_uw[4], dS_c[4], dS_dw[4];
  double  A[4][4], Ma[4][4], Mb[4][4], Mc[4][4], V0[4];
  double  M_dw[4][4];
  
  Stencil  *st;
  LongChar *lc = NULL;

  /* --- Boundary conditions in row at the top or bottom -- --------- */

  if (to_observer) {
    row_left = (geometry->Nz-1) * geometry->Nx;
    bvalue = geometry->bvalue[BOTTOM];
  } else {
    row_left = 0;
    bvalue = geometry->bvalue[TOP];
  }
  row_right = row_left + geometry->Nx;

  switch (bvalue) {
  case IRRADIATED:
    if (to_observer) {
      for (l = 0;  l < geometry->Nx;  l++)
	I[0][row_left + l] = geometry->Ibottom[nspect][l];
    } else {
      for (l = 0;  l < geometry->Nx;  l++)
	I[0][row_left + l] = geometry->Itop[nspect][l];
    }
    break;
  case ZERO:
    for (l = row_left;  l < row_right;  l++) I[0][l] = 0.0;
    break;
  case THERMALIZED:
    if (to_observer) {
      local = (geometry->Nz-1) * geometry->Nx;
      for (l = 0;  l < geometry->Nx;  l++) {
	st = &geometry->stencil[mu][local];

        T_uw = Quadr(atmos.T, st, UPWIND, monotonic);    
	Planck(1, &T_uw, spectrum.lambda[nspect], &Bnu_uw);
	Planck(1, &atmos.T[local], spectrum.lambda[nspect], &Bnu);

	chi_uw  = Quadr(chi, st, UPWIND, monotonic);
	dtau_uw = 0.5 * (chi_uw + chi[local]) * st->ds[UPWIND];

	I[0][local] = Bnu - (Bnu_uw - Bnu) / dtau_uw;
        local++;
      }
    } else {
      Error(ERROR_LEVEL_2, routineName, 
	    "Boundary condition THERMALIZED not implemented for TOP");
    }
    break;

  }
  if (Psi)
    for (l = row_left;  l < row_right;  l++) Psi[l] = 0.0;

  /* --- Assume irradiation is unpolarized --          -------------- */

  for (n = 1;  n < 4;  n++) {
    for (l = row_left;  l < row_right;  l++) I[n][l] = 0.0;
  } 

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
  /* --- Finally, go through the grid to solve the non-local transfer
         equation --                                   -------------- */

  for (k = kstart;  k != kend + dk;  k += dk) {
    lc = &geometry->longchar[mu][sweep*geometry->Nz + k];
    lstart = lc->lstart;
 
    for (lp = 0, l = lstart;  lp < geometry->Nx;  lp++, l += dl) {
      local = k*geometry->Nx + MODULO(l, geometry->Nx);
      st    = &geometry->stencil[mu][local];

      if (l == lc->lstart  &&  lc->Nlc > 0)
	SolveLongStokes(nspect, lc, local, chi, S, I, I_uw);
      else {
	for (n = 0;  n < 4;  n++)
	  I_uw[n] = Quadr(I[n], st, upwind, monotonic);
      }

      chi_uw = Quadr(chi, st, upwind, monotonic);
      dchi_uw = (chi[local] - chi_uw) / st->ds[upwind];
      StokesK_2D(nspect, st, upwind, monotonic, chi_uw, K_uw);
      for (n = 0;  n < 4;  n++) {
	S_uw[n] = Quadr(S[n], st, upwind, monotonic);
      }
      StokesK(nspect, local, chi[local], K_loc);
      	
      if (k == kend) {

	dtau_uw = 0.5 * (chi[local] + chi_uw) * st->ds[upwind];
	
	/* --- Piecewise linear at end of ray --       -------------- */

	for (n = 0;  n < 4;  n++) {
	  dS_uw[n] = (S_loc[n] - S_uw[n]) / dtau_uw;
	  dS_c[n]  = dS_uw[n];
	
	  for(m = 0;  m < 4;  m++) {
	    dK_uw[n][m] = (K_loc[n][m] - K_uw[n][m]) / dtau_uw;
	    dK_c[n][m]  = dK_uw[n][m];
	  }
	}
      } else {

	/* --- Bezier piecewise elsewhere --          --------------- */

	chi_dw  = Quadr(chi, st, downwind, monotonic);
	dchi_dw = (chi_dw - chi[local]) / st->ds[downwind];
	StokesK_2D(nspect, st, downwind, monotonic, chi_dw, K_dw);
	
	dchi_c  = cent_deriv(st->ds[upwind], st->ds[downwind],
			     chi_uw, chi[local], chi_dw);

	/* --- Upwind optical path length --           -------------- */
	
	ds03 = st->ds[upwind] * 0.3333333333333333333333333;
	c1   = chi_uw     + ds03 * dchi_uw;
	c2   = chi[local] - ds03 * dchi_c;
	dtau_uw = st->ds[upwind] * (chi[local] + chi_uw +
				    c1 + c2) * 0.25;
      
	/* --- Downwind optical path length --          ------------- */

	ds03 = st->ds[downwind] * 0.3333333333333333333333333;
	c1   = chi_dw     - ds03 * dchi_dw;
	c2   = chi[local] + ds03 * dchi_c;
	dtau_dw = st->ds[downwind] * (chi[local] + chi_dw +
				      c1 + c2) * 0.25;
    
	for (n = 0;  n < 4;  n++) {
	  S_dw[n]  = Quadr(S[n], st, downwind, monotonic);
	  S_loc[n] = S[n][local];
	  dS_uw[n] = (S_loc[n] - S_uw[n]) / dtau_uw;

	  for(m = 0;  m < 4;  m++)
	    dK_uw[n][m] = (K_loc[n][m] - K_uw[n][m]) / dtau_uw;
	}

	cent_deriv_vec(dS_c, dtau_uw, dtau_dw, S_uw, S_loc, S_dw);
	cent_deriv_mat(dK_c, dtau_uw, dtau_dw, K_uw, K_loc, K_dw);
      }
      m4m(K_uw, K_uw, Ma);  // Ku # Ku
      m4m(K_loc, K_loc, A); // K0 # K0

      /* --- Compute interpolation parameters --     ---------------- */

      Bezier3_coeffs(dtau_uw, &alpha, &beta, &gamma, &theta, &eps);

      /* --- Diagonal operator --                    ---------------- */

      if (Psi) Psi[local] = beta + theta;

      dt03 = dtau_uw * 0.3333333333333333333333333;
      for(j = 0;  j < 4;  j++){
	for(i = 0;  i < 4;  i++){
	  M_dw[j][i] = ident[j][i] + beta * K_loc[j][i] + theta *
	    (dt03 * (A[j][i] + dK_c[j][i] + K_loc[j][i]) + K_loc[j][i]);
	    
	  Ma[j][i] = eps * ident[j][i] - alpha * K_uw[j][i] + gamma *
	    (dt03 * (Ma[j][i] - dK_uw[j][i] + K_uw[j][i]) - K_uw[j][i]);
	  
	  Mb[j][i] = alpha * ident[j][i] + gamma * (ident[j][i] -
						    dt03 * K_uw[j][i]);
	  Mc[j][i] = beta* ident[j][i] + theta * (ident[j][i] +
						   dt03 * K_loc[j][i]);
	}
      }
      
      /* --- Here I am doing Ma*stk + Mb * Su + Mc * S0 + 
	 (gam * dS0 - theta * dSu) * dtau / 3.0 to compute the 
	 right-hand term
           --                                      ------------------ */

      memset(V0, 0, 4*sizeof(double));
	
      for(i = 0;  i < 4;  i++){
	for(j = 0;  j < 4;  j++){
	  V0[i] += Ma[i][j] * I_uw[j] + Mb[i][j] * S_uw[j] +
	    Mc[i][j] * S_loc[j];
	}
	V0[i] += dt03 * (gamma * dS_c[i] - theta * dS_uw[i]);
      }
      /* --- Solve linear system to get the intensity -- ----------- */
      
      solveLinearFast(M_dw, V0);

      /* --- Finally, store intensities --           --------------- */
      
      for(n = 0;  n < 4;  n++) I[n][local] = V0[n];
    }
  }
}
/* ------- end ------------------------- Piece_Stokes_Bezier3_2D.c -- */


/* ------- begin -------------------------- Piecewise_Bezier3_2D.c -- */

void Piecewise_Bezier3_2D(Geometry *geometry, int nspect, int mu,
			  bool_t to_observer, double *chi, double *S,
			  double *I, double *Psi)
{

  /* --- Cubic short-char Bezier solver.
         Coded by J. de la Cruz Rodriguez (ISP-SU 2017).
         Adapted to 2-D version by Han Uitenbroek.

         Reference(s):
         J. de la Cruz Rodriguez & N. Piskunov (2013)
         --                                        ------------------ */
  
  const char routineName[] = "Piecewise_Bezier3_2D";
  register int l, k, n;

  bool_t  monotonic = TRUE;
  enum    boundval bvalue;
  enum    direction upwind, downwind;
  enum    sweep sweep;
  int     Nx = geometry->Nx, Nz = geometry->Nz, local, dl, lstart;
  int     Nrow, row_left, row_right, dk, kstart, kend;
  double  mux = geometry->mux[mu], I_uw, *Iboundary, dtau_uw, dtau_dw;
  double  c1, c2, dS_uw, dS_dw, S_uw, chi_uw, S_dw, chi_dw, dt03, w[2];
  double  eps = 0, alpha = 0, beta = 0, gamma = 0, theta = 0;
  double  dchi_dw, dchi_uw, dchi_c, dS_c, ds03;
  double  Bnu_uw, Bnu, T_uw;
  
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

      if (geometry->hboundary == PERIODIC &&
	  l == lc->lstart && lc->Nlc > 0)
	I_uw = SolveLong(lc, local, chi, S, I);
      else
	I_uw = Quadr(I, st, upwind, monotonic);
      
      chi_uw  = Quadr(chi, st, upwind, monotonic);
      dchi_uw = (chi[local] - chi_uw) / st->ds[upwind];
      S_uw    = Quadr(S, st, upwind, monotonic);

      if (k == kend ||
	  (geometry->hboundary == FIXED && n == Nrow-1)) {
	
	dtau_uw = 0.5 * (chi[local] + chi_uw) * st->ds[upwind];
       
	/* --- Piecewise linear at end of ray --       -------------- */

	dS_c = (S[local] - S_uw) / dtau_uw;
      } else {

	/* --- Bezier integration --                   -------------- */
	
	chi_dw  = Quadr(chi, st, downwind, monotonic);
	dchi_dw = (chi_dw - chi[local]) / st->ds[downwind];
	
	dchi_c  = cent_deriv(st->ds[upwind], st->ds[downwind],
			     chi_uw, chi[local], chi_dw);
	
	/* --- Upwind optical path length --          --------------- */

	ds03 = st->ds[upwind] * 0.333333333333333333333333333;
	c1   = chi[local] - ds03 * dchi_c;
	c2   = chi_uw     + ds03 * dchi_uw;
	dtau_uw = st->ds[upwind] * (chi[local] + chi_uw +
				    c1 + c2) * 0.25;

 	/* --- Downwind optical path length --          ------------- */

	ds03 = st->ds[downwind] * 0.333333333333333333333333333;
	c1   = chi[local] + ds03 * dchi_c;
	c2   = chi_dw     - ds03 * dchi_dw;
	dtau_dw = st->ds[downwind] * (chi[local] + chi_dw +
				      c1 + c2) * 0.25;
    
	/* --- dS/dt at central point --               -------------- */

	S_dw = Quadr(S, st, downwind, monotonic);
	dS_c = cent_deriv(dtau_uw, dtau_dw, S_uw, S[local], S_dw);
      }
      /* --- Source function control points --       ---------------- */

      dt03  = dtau_uw * 0.333333333333333333333333333;
      dS_uw = (S[local] - S_uw) / dtau_uw;
      c1    = S[local] - dt03 * dS_c;
      c2    = S_uw     + dt03 * dS_uw;       
     
      /* --- Compute interpolation parameters --     ---------------- */
      
      Bezier3_coeffs(dtau_uw, &alpha, &beta, &gamma, &theta, &eps);

      /* --- Solve integral in this interval --      ---------------- */
       
      I[local] = I_uw*eps + beta*S[local] + alpha*S_uw +
	theta * c1 + gamma * c2; 

      /* --- Diagonal operator --                    ---------------- */

      if (Psi) Psi[local] = beta + theta;
    }
  }
}
/* ------- end ---------------------------- Piecewise_Bezier3_2D.c -- */
