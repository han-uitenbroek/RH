/* ------- file: -------------------------- bezier_3D.c -- ----------

       Version:       rh2.0, 3-D Cartesian, short characteristics

   Cubic DELO-Bezier (polarized) and cubic short-char Bezier solvers.
   
   References: de la Cruz Rodriguez & Piskunov (2013), Auer (2003)
               (Derivatives) Fritsch & Butland (1984),
	       
   1D version Coded by J. de la Cruz Rodriguez (ISP-SU 2017),
   adapted to 3D by Han Uitenbroek

       Last modified: Wed Feb 20 13:52:27 2019 --

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
  S_uw    = Interpolate_3D(S, geometry, st_uw, l, m);

  if (k == kend) {
    dtau_uw = 0.5 * (chi_uw + chi[local]) * st_uw->ds;
    dS_uw   = (S_uw - S[local]) / dtau_uw;
    w2(dtau_uw, w);

    /* --- Piecewise linear integration in last layer -- ------------ */

    I[local] = I_uw*(1.0 - w[0]) + w[0]*S[local] + w[1]*dS_uw;
    if (Psi) Psi[local] = w[0] - w[1]/dtau_uw;
    
  } else {


    /* --- The downwind quantities --                  -------------- */

    chi_dw  = Interpolate_3D(chi, geometry, st_dw, l, m);
    dchi_uw =  (chi[local] - chi_uw) / st_uw->ds;
    dchi_dw =  (chi_dw - chi[local]) / st_dw->ds;
    dchi_c  = cent_deriv(st_uw->ds, st_dw->ds,
			 chi_uw, chi[local], chi_dw);

    /* --- Upwind optical path length --               ------------- */

    ds03 = st_uw->ds / 3.0;
    c1   = MAX(chi[local] - ds03 * dchi_c,  0.0);
    c2   = MAX(chi_uw     + ds03 * dchi_uw, 0.0);
    dtau_uw = st_uw->ds * (chi[local] + chi_uw + c1 + c2) * 0.25;

    /* --- Downwind optical path length --             ------------- */

    ds03 = st_dw->ds / 3.0;
    c1   = MAX(chi[local] + ds03 * dchi_c,  0.0);
    c2   = MAX(chi_dw     - ds03 * dchi_dw, 0.0);
    dtau_dw = st_dw->ds * (chi[local] + chi_dw + c1 + c2) * 0.25;

    /* --- dS/dt at central point --                  -------------- */

    S_dw = Interpolate_3D(S, geometry, st_dw, l, m);
    dS_c = cent_deriv(dtau_uw, dtau_dw, S_uw, S[local], S_dw);

    /* --- Source function control points --          -------------- */
       
    dS_uw = (S[local] - S_uw) / dtau_uw;
    
    dt03  = dtau_uw / 3.0;
    c1    = MAX(S[local] - dt03 * dS_c , 0.0);
    c2    = MAX(S_uw     + dt03 * dS_uw, 0.0);       
     
    /* --- Compute interpolation parameters --        -------------- */
       
    Bezier3_coeffs(dtau_uw, &beta, &alpha, &theta, &gamma, &eps);

    /* --- Solve integral in this interval --         -------------- */
       
    I[local] = I_uw*eps + alpha*S[local] + beta*S_uw +
      gamma * c1 + theta * c2; 

    /* --- Diagonal operator --                       -------------- */

    if (Psi) Psi[local] = alpha + gamma;
  }
}
/* ------- end ---------------------------- Piecewise_Bezier3_3D.c -- */

/* ------- begin -------------------------- Piece_Stokes_Stokes_3D.c  */

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
  double dK_uw[4][4], dK_c[4][4], dK_dw[4][4];
  double A[4][4], Ma[4][4], Mb[4][4], Mc[4][4], V0[4], V1[4];
  float  M_dw[4][4];

  local = k*geometry->Nplane + m*geometry->Nx + l;

  /* --- The upwind quantities --                      -------------- */

  chi_uw = Interpolate_3D(chi, geometry, st_uw, l, m);
  StokesK_3D(nspect, geometry, st_uw, l, m, chi_uw, K_uw);
  StokesK(nspect, local, chi[local], K_loc);

  for (n = 0;  n < 4;  n++) {
    S_uw[n]  = Interpolate_3D(S[n], geometry, st_uw, l, m);
  }

  if (k == kend) {
	
    /* --- Piecewise linear at end of ray --           -------------- */

    dtau_uw = 0.5 * (chi_uw + chi[local]) * st_uw->ds;
    for (n = 0;  n < 4;  n++) {
      dS_uw[n] = (S[n][local] - S_uw[n]) / dtau_uw;
    }
    w2(dtau_uw, w);

    for (n = 0;  n < 4;  n++) {
      V0[n] = w[0]*S[n][local] + w[1]*dS_uw[n];
      for (p = 0;  p < 4;  p++) {
	A[n][p]    = -w[1]/dtau_uw * K_uw[n][p];
	M_dw[n][p] = (w[0] - w[1]/dtau_uw) * K_loc[n][p];
      }
      A[n][n]    = 1.0 - w[0];
      M_dw[n][n] = 1.0;
    }
    for (n = 0;  n < 4;  n++)
      for (p = 0;  p < 4;  p++) 
	V0[n] += A[n][p] * I_uw[p];

    /* --- Solve linear system --                      -------------- */
  
    SIMD_MatInv(M_dw[0]);   // Invert Md
    m4v(M_dw, V0, V1);      // Multiply Md^-1 * V0
  
    for (n = 0;  n < 4;  n++) I[n][local] = V1[n];
	
    if (Psi) Psi[local] = w[0] - w[1] / dtau_uw;
   } else {

    /* --- Bezier piecewise elsewhere --               -------------- */

    chi_dw  = Interpolate_3D(chi, geometry, st_dw, l, m);
    dchi_uw = (chi[local] - chi_uw) / st_uw->ds;
    dchi_dw = (chi_dw - chi[local]) / st_dw->ds;
    dchi_c  = cent_deriv(st_uw->ds, st_dw->ds,
			 chi_uw, chi[local], chi_dw);

    /* --- Upwind optical path length --               -------------- */

    ds03 = st_uw->ds / 3.0;
    c1   = MAX(chi[local] - ds03 * dchi_c,  0.0);
    c2   = MAX(chi_uw     + ds03 * dchi_uw, 0.0);
    dtau_uw = st_uw->ds * (chi[local] + chi_uw + c1 + c2) * 0.25;

    for (n = 0;  n < 4;  n++) {
      for(p = 0;  p < 4;  p++)
	dK_uw[n][p] = (K_loc[n][p] - K_uw[n][p]) / dtau_uw;
    }

    /* --- Downwind optical path length --             ------------- */

    ds03 = st_dw->ds / 3.0;
    c1   = MAX(chi[local] + ds03 * dchi_c,  0.0);
    c2   = MAX(chi_dw     - ds03 * dchi_dw, 0.0);
    dtau_dw = st_dw->ds * (chi[local] + chi_dw + c1 + c2) * 0.25;
    
    StokesK_3D(nspect, geometry, st_dw, l, m, chi_dw, K_dw);

    for (n = 0;  n < 4;  n++) {
      S_dw[n]  = Interpolate_3D(S[n], geometry, st_dw, l, m);
      S_loc[n] = S[n][local];
      dS_uw[n] = (S_loc[n] - S_uw[n]) / dtau_uw;

      for(p = 0;  p < 4;  p++)
	dK_dw[n][p] = (K_dw[n][p] - K_loc[n][p]) / dtau_dw;
    }

    cent_deriv_mat(dK_c, dtau_uw, dtau_dw, K_uw, K_loc, K_dw);
    cent_deriv_vec(dS_c, dtau_uw, dtau_dw, S_uw, S_loc, S_dw);

    m4m(K_uw, K_uw, Ma);  // Ku # Ku
    m4m(K_loc, K_loc, A); // K0 # K0

    /* --- Compute interpolation parameters --         -------------- */
       
    Bezier3_coeffs(dtau_uw, &beta, &alpha, &theta, &gamma, &eps);

    dt03 = dtau_uw / 3.0;
    for(j = 0;  j < 4;  j++){
      for(i = 0;  i < 4;  i++){
	M_dw[j][i] = ident[j][i] + alpha * K_loc[j][i] - gamma *
	  -(dt03 * (A[j][i] + dK_c[j][i] + K_loc[j][i]) + K_loc[j][i]);
	
	Ma[j][i] = eps * ident[j][i] - beta * K_uw[j][i] + theta *
	  (dt03 * (Ma[j][i] + dK_uw[j][i] + K_uw[j][i]) - K_uw[j][i]);
	  
	Mb[j][i] = beta * ident[j][i] + theta * (ident[j][i] -
						 dt03 * K_uw[j][i]);
	Mc[j][i] = alpha* ident[j][i] + gamma * (ident[j][i] +
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
    /* --- Solve linear system to get the intensity -- -------------- */
      
    SIMD_MatInv(M_dw[0]);   // Invert Md
    m4v(M_dw, V0, V1);      // Multiply Md^-1 * V0

    /* --- Finally, store intensities --               -------------- */
	
    for(n = 0;  n < 4;  n++) I[n][local] = V1[n];
	
    /* --- Diagonal operator --                        -------------- */

    if (Psi) Psi[local] = alpha + gamma;
  }
}
/* ------- end ---------------------------- Piece_Stokes_Bezier3_3D.c */
