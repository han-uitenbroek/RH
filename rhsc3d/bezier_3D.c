/* ------- file: -------------------------- bezier_3D.c -- ----------

       Version:       rh2.0, 3-D Cartesian, short characteristics

   Cubic DELO-Bezier (polarized) and cubic short-char Bezier solvers.
   
   References: de la Cruz Rodriguez & Piskunov (2013), Auer (2003)
               (Derivatives) Fritsch & Butland (1984),
	       
   1D version Coded by J. de la Cruz Rodriguez (ISP-SU 2017),
   adapted to 3D by Han Uitenbroek

       Last modified: Mon Jan  6 14:59:35 2020 --

       --------------------------                      ----------RH-- */

#include <math.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "geometry.h"
#include "spectrum.h"
#include "inputs.h"
#include "error.h"
#include "bezier.h"


/* --- Function prototypes --                          -------------- */


/* --- Identity matrix --                              -------------- */

static const double ident[4][4] =
  {{1.0, 0.0, 0.0, 0.0},
   {0.0, 1.0, 0.0, 0.0},
   {0.0, 0.0, 1.0, 0.0},
   {0.0, 0.0, 0.0, 1.0}};


/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern Spectrum spectrum;
extern InputData input;
extern char messageStr[];

/* ------- begin -------------------------- Piece_Stokes_Bezier3_3D.c */

void Piece_Stokes_Bezier3_3D(Geometry *geometry, Stencil *st_uw,
			     Stencil *st_dw,
			     int nspect, int k, int kend, int l, int m,
			     double *I_uw,
			     double *chi, double **S, double **I, double *Psi)
{
  /* --- Cubic short-char Bezier Polarization solver.
         Coded by J. de la Cruz Rodriguez (ISP-SU 2017).
         Adapted to 3-D version by Han Uitenbroek.

         Reference(s):
         J. de la Cruz Rodriguez & N. Piskunov (2013)
         --                                        ------------------ */

  register int n, p, i, j;

  int    local;
  double chi_uw, chi_dw, S_uw[4], S_dw[4], dS_dw[4], dS_uw[4], w[2];
  double eps = 0, alpha = 0, beta = 0, gamma = 0, theta = 0;
  double dtau_uw, dtau_dw, S_loc[4], dS_c[4], dchi_uw, dchi_dw, dchi_c;
  double K_loc[4][4], K_uw[4][4], K_dw[4][4], c1, c2, ds03, dt03;
  double dK_uw[4][4], dK_c[4][4];
  double A[4][4], Ma[4][4], Mb[4][4], Mc[4][4], V0[4];
  double M_dw[4][4];

  local = k*geometry->Nplane + m*geometry->Nx + l;

  /* --- The upwind quantities --                      -------------- */
  
  chi_uw  = Interpolate_3D(chi, geometry, st_uw, l, m);
  dchi_uw = (chi[local] - chi_uw) / st_uw->ds;
  StokesK_3D(nspect, geometry, st_uw, l, m, chi_uw, K_uw);
  for (n = 0;  n < 4;  n++) {
    S_uw[n]  = Interpolate_3D(S[n], geometry, st_uw, l, m);
  }
  StokesK(nspect, local, chi[local], K_loc);
  
  if (k == kend) {

    /* --- Piecewise linear at end of ray --           -------------- */

    dtau_uw = 0.5 * (chi[local] + chi_uw) * st_uw->ds;
    
    for (n = 0;  n < 4;  n++) {
      dS_uw[n] = (S_loc[n] - S_uw[n]) / dtau_uw;
      dS_c[n]  = dS_uw[n];
      
      for(p = 0;  p < 4;  p++) {
	dK_uw[n][p] = (K_loc[n][p] - K_uw[n][p]) / dtau_uw;
	dK_c[n][p]  = dK_uw[n][p];
      }
    }
  } else {

    /* --- Bezier piecewise elsewhere --               -------------- */

    chi_dw  = Interpolate_3D(chi, geometry, st_dw, l, m);
    dchi_dw = (chi_dw - chi[local]) / st_dw->ds;
    StokesK_3D(nspect, geometry, st_dw, l, m, chi_dw, K_dw);

    dchi_c  = cent_deriv(st_uw->ds, st_dw->ds,
			 chi_uw, chi[local], chi_dw);

    /* --- Upwind optical path length --               -------------- */
  
    ds03 = st_uw->ds * 0.3333333333333333333333333;
    c1   = chi_uw     + ds03 * dchi_uw;
    c2   = chi[local] - ds03 * dchi_c;
    dtau_uw = st_uw->ds * (chi[local] + chi_uw + c1 + c2) * 0.25;

    /* --- Downwind optical path length --             ------------- */

    ds03 = st_dw->ds * 0.3333333333333333333333333;
    c1   = chi_dw     - ds03 * dchi_dw;
    c2   = chi[local] + ds03 * dchi_c;
    dtau_dw = st_dw->ds * (chi[local] + chi_dw + c1 + c2) * 0.25;
    
    for (n = 0;  n < 4;  n++) {
      S_dw[n]  = Interpolate_3D(S[n], geometry, st_dw, l, m);
      S_loc[n] = S[n][local];
      dS_uw[n] = (S_loc[n] - S_uw[n]) / dtau_uw;

      for(p = 0;  p < 4;  p++)
	dK_uw[n][p] = (K_loc[n][p] - K_uw[n][p]) / dtau_uw;
    }

    cent_deriv_vec(dS_c, dtau_uw, dtau_dw, S_uw, S_loc, S_dw);
    cent_deriv_mat(dK_c, dtau_uw, dtau_dw, K_uw, K_loc, K_dw);
  }
  m4m(K_uw, K_uw, Ma);  // Ku # Ku
  m4m(K_loc, K_loc, A); // K0 # K0

  /* --- Compute interpolation parameters --         ---------------- */
       
  Bezier3_coeffs(dtau_uw, &alpha, &beta, &gamma, &theta, &eps);

  /* --- Diagonal operator --                          -------------- */
  
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
           --                                         --------------- */
    
  memset(V0, 0, 4*sizeof(double));
  for(i = 0;  i < 4;  i++){
    for(j = 0;  j < 4;  j++){
      V0[i] += Ma[i][j] * I_uw[j] + Mb[i][j] * S_uw[j] +
	Mc[i][j] * S_loc[j];
    }
    V0[i] += dt03 * (gamma * dS_c[i] - theta * dS_uw[i]);
  }
  /* --- Solve linear system to get the intensity -- ---------------- */
      
  solveLinearFast(M_dw, V0);

  /* --- Finally, store intensities --           -------------------- */
    
  for(n = 0;  n < 4;  n++) I[n][local] = V0[n];
}
/* ------- end ---------------------------- Piece_Stokes_Bezier3_3D.c */

/* ------- begin -------------------------- Piecewise_Bezier3_3D.c -- */

void Piecewise_Bezier3_3D(Geometry *geometry, Stencil *st_uw, Stencil *st_dw,
			  int k, int kend, int l, int m, double I_uw,
			  double *chi, double *S, double *I, double *Psi)
{
  /* --- Cubic short-char Bezier solver.
         Coded by J. de la Cruz Rodriguez (ISP-SU 2017).
         Adapted to 3-D version by Han Uitenbroek.

         Reference(s):
         J. de la Cruz Rodriguez & N. Piskunov (2013)
         --                                        ------------------ */
  
  int    local;
  double  dtau_uw, dtau_dw;
  double  c1, c2, dS_uw, dS_dw, S_uw, chi_uw, S_dw, chi_dw, dt03, w[2];
  double  eps = 0, alpha = 0, beta = 0, gamma = 0, theta = 0;
  double  dchi_dw, dchi_uw, dchi_c, dS_c, ds03;

  local = k*geometry->Nplane + m*geometry->Nx + l;

  /* --- The upwind quantities --                      -------------- */

  chi_uw  = Interpolate_3D(chi, geometry, st_uw, l, m);
  dchi_uw =  (chi[local] - chi_uw) / st_uw->ds;
  S_uw    = Interpolate_3D(S, geometry, st_uw, l, m);

  if (k == kend) {

    dtau_uw = 0.5 * st_uw->ds * (chi[local] + chi_uw);
    
    /* --- Piecewise linear at end of ray --           -------------- */

    dS_c = (S[local] - S_uw) / dtau_uw;
   } else {

    /* --- Bezier interpolation --                     -------------- */
    
    chi_dw  = Interpolate_3D(chi, geometry, st_dw, l, m);
    dchi_dw = (chi_dw - chi[local]) / st_dw->ds;

    dchi_c  = cent_deriv(st_uw->ds, st_dw->ds,
			 chi_uw, chi[local], chi_dw);

    /* --- Upwind optical path length --             ---------------- */

    ds03 = st_uw->ds * 0.333333333333333333333333333;
    c1   = chi[local] - ds03 * dchi_c;
    c2   = chi_uw     + ds03 * dchi_uw;
    dtau_uw = st_uw->ds * (chi[local] + chi_uw + c1 + c2) * 0.25;

    /* --- Downwind optical path length --             ------------- */

    ds03 = st_dw->ds * 0.333333333333333333333333333;
    c1   = chi[local] + ds03 * dchi_c;
    c2   = chi_dw     - ds03 * dchi_dw;
    dtau_dw = st_dw->ds * (chi[local] + chi_dw + c1 + c2) * 0.25;
    S_dw   = Interpolate_3D(S, geometry, st_dw, l, m);

    /* --- dS/dt at central point --                   -------------- */

    dS_c = cent_deriv(dtau_uw, dtau_dw, S_uw, S[local], S_dw);
  }
  /* --- Source function control points --             -------------- */
       
  dt03  = dtau_uw * 0.333333333333333333333333333;
  dS_uw = (S[local] - S_uw) / dtau_uw;
  c1    = S[local] - dt03 * dS_c;
  c2    = S_uw     + dt03 * dS_uw;       
     
  /* --- Compute interpolation parameters --           -------------- */
       
  Bezier3_coeffs(dtau_uw, &alpha, &beta, &gamma, &theta, &eps);
    
  /* --- Solve integral in this interval --            -------------- */
  
  I[local] = I_uw*eps + beta*S[local] + alpha*S_uw +
    theta * c1 + gamma * c2; 
    
  /* --- Diagonal operator --                          -------------- */

    if (Psi) Psi[local] = beta + theta;
}
/* ------- end ---------------------------- Piecewise_Bezier3_3D.c -- */
