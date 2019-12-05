/* ------- file: -------------------------- bezier_1D.c -------------

   Cubic DELO-Bezier (polarized) and cubic short-char Bezier solvers.
   
   References: de la Cruz Rodriguez & Piskunov (2013), Auer (2003)
               (Derivatives) Fritsch & Butland (1984),
	       
   Coded by J. de la Cruz Rodriguez (ISP-SU 2017)

   Modifications:
           2017-03-12, JdlCR: Created!

       Last modified: Wed Feb 20 13:30:27 2019 --

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


/* ------- begin -------------------------- Piece_Stokes_Bezier3_1D.c */

void Piece_Stokes_Bezier3_1D(int nspect, int mu, bool_t to_obs,
			     double *chi, double **S, double **I,
			     double *Psi)
{
  /* --- Cubic DELO-Bezier solver for polarized light
         Coded by J. de la Cruz Rodriguez (ISP-SU 2017)

         Reference(s):
         J. de la Cruz Rodriguez & N. Piskunov (2013)
         --                                        ------------------ */
  
  const char routineName[] = "Piece_Stokes_Bezier3_1D";
  register int k, n, m, i, j;
  
  int    Ndep = geometry.Ndep, k_start, k_end, dk;
  double dtau_uw, dtau_dw = 0.0, c1, c2, w[3], dsdn2, dchi_dn,
         I_upw[4], Bnu[2];
  double dchi_up,dchi_c,dt03;
  double dsup,dsdn,dt,eps=0,alpha=0,beta=0,gamma=0,theta=0;
  double Ku[4][4], K0[4][4], Kd[4][4], dKu[4][4], dK0[4][4];
  double Su[4], S0[4], Sd[4], dSu[4], dS0[4];
  double A[4][4], Ma[4][4], Mb[4][4], Mc[4][4], V0[4], V1[4];
  double imu = 1.0 / geometry.muz[mu];
  float Md[4][4];
  double *z = geometry.height;

  
  if (to_obs) {
    dk      = -1;
    k_start = Ndep-1;
    k_end   = 0;
  } else {
    dk      = 1;
    k_start = 0;
    k_end   = Ndep-1;
  }
  dtau_uw = 0.5 * imu * (chi[k_start] + chi[k_start+dk]) *
    fabs(z[k_start] - z[k_start+dk]);
  
  /* --- Boundary conditions --                        -------------- */

  if (to_obs) {
    switch (geometry.vboundary[BOTTOM]) {
    case ZERO:
      for (n = 0;  n < 4;  n++) I_upw[n] = 0.0;
      break;
    case THERMALIZED:
      Planck(2, &atmos.T[Ndep-2], spectrum.lambda[nspect], Bnu);
      I_upw[0] = Bnu[1] - (Bnu[0] - Bnu[1]) / dtau_uw;
      for (n = 1;  n < 4;  n++) I_upw[n] = 0.0;
      break;
    case IRRADIATED:
      I_upw[0] = geometry.Ibottom[nspect][mu];
      for (n = 1;  n < 4;  n++) I_upw[n] = 0.0;
      break;
    case REFLECTIVE:
      sprintf(messageStr, "Boundary condition not implemented: %d",
	      geometry.vboundary[BOTTOM]);
      Error(ERROR_LEVEL_2, routineName, messageStr);
    }
  } else {
    switch (geometry.vboundary[TOP]) {
    case ZERO:
      for (n = 0;  n < 4;  n++) I_upw[n] = 0.0;
      break;
    case IRRADIATED:
      I_upw[0] = geometry.Itop[nspect][mu];
      for (n = 1;  n < 4;  n++) I_upw[n] = 0.0;
      break;
    default:
      sprintf(messageStr, "Boundary condition not implemented: %d",
	      geometry.vboundary[TOP]);
      Error(ERROR_LEVEL_2, routineName, messageStr);
    }
  }

  for (n = 0;  n < 4;  n++) I[n][k_start] = I_upw[n];
  if (Psi) Psi[k_start] = 0.0;
  
  k=k_start+dk;
  dsup = fabs(z[k] - z[k-dk]) * imu;
  dsdn = fabs(z[k+dk] - z[k]) * imu;
  dchi_up= (chi[k] - chi[k-dk])/dsup;

  
  /* ---  dchi/ds at central point--               ------------------ */
  
  dchi_c = cent_deriv(dsup,dsdn,chi[k-dk],chi[k],chi[k+dk]);
  
  
  /* --- Upwind path_length (BEzier3 integration) -- ---------------- */

  c2 = MAX(chi[k]    - (dsup/3.0) * dchi_c,  0.0);
  c1 = MAX(chi[k-dk] + (dsup/3.0) * dchi_up, 0.0);
  
  dtau_uw = 0.25 * dsup * (chi[k] + chi[k-dk] + c1 + c2);

  
  /* --- Ku, K0 and dKu, dSu -                     ------------------ */
  
  StokesK(nspect, k_start,    chi[k_start],    Ku);
  StokesK(nspect, k_start+dk, chi[k_start+dk], K0);

  Svec(k_start,    S, Su);
  Svec(k_start+dk, S, S0);

  /* --- Assume side derivative in the first interval -- ------------ */
  
  for(n = 0;  n < 4;  n++){
    dSu[n] = (S0[n] - Su[n]) / dtau_uw;
    
    for(m = 0;  m < 4;  m++)
      dKu[n][m] = (K0[n][m] - Ku[n][m]) / dtau_uw;
  }
  
  /* --- Solve transfer along ray --                   -------------- */

  for (k = k_start+dk;  k != k_end;  k += dk) {      
      
    /* --- dchi/ds at downwind point --                -------------- */
      
    dsdn = fabs(z[k+dk] - z[k]) * imu;
      
    if(fabs(k - k_end) > 1){
      dsdn2 = fabs(z[k+2*dk] - z[k+dk]) * imu;
      dchi_dn = cent_deriv(dsdn, dsdn2, chi[k], chi[k+dk], chi[k+2*dk]);       
    } else
      dchi_dn = (chi[k+dk] - chi[k])/dsdn;      
      
    /* --- Make sure that c1 and c2 don't do below zero -- ---------- */
      
    c2 = MAX(chi[k]    + (dsdn/3.0) * dchi_c , 0.0);
    c1 = MAX(chi[k+dk] - (dsdn/3.0) * dchi_dn, 0.0);
          
    /* --- Bezier3 integrated dtau --              ------------------ */
      
    dtau_dw = 0.25 * dsdn * (chi[k] + chi[k+dk] + c1 + c2);
    dt = dtau_uw, dt03 = dt / 3.0;
  
    /* --- Bezier3 coeffs. --                      ------------------ */
      
    Bezier3_coeffs(dt, &beta, &alpha, &theta, &gamma, &eps);
   
    /* --- Diagonal operator --                    ------------------ */
      
    if(Psi) Psi[k] = alpha + gamma;
   
    /* ---- get algebra in place --                ------------------ */
      
    StokesK(nspect, k+dk, chi[k+dk], Kd);
    Svec(k+dk, S, Sd);

    cent_deriv_mat(dK0, dtau_uw, dtau_dw, Ku, K0, Kd);
    cent_deriv_vec(dS0, dtau_uw, dtau_dw, Su, S0, Sd);

    m4m(Ku, Ku, Ma); // Ku # Ku
    m4m(K0, K0, A ); // K0 # K0

    for(j = 0;  j < 4;  j++){
      for(i = 0;  i < 4;  i++){
	Md[j][i] = ident[j][i] + alpha * K0[j][i] - gamma *
	  -(dt03 * (A[j][i] + dK0[j][i] + K0[j][i]) + K0[j][i]);
	  
	Ma[j][i] = eps * ident[j][i] - beta * Ku[j][i] + theta *
	  (dt03 * (Ma[j][i] + dKu[j][i] + Ku[j][i]) - Ku[j][i]);
	  
	Mb[j][i] = beta * ident[j][i] + theta * (ident[j][i] -
						 dt03 * Ku[j][i]);
	Mc[j][i] = alpha* ident[j][i] + gamma * (ident[j][i] +
						 dt03 * K0[j][i]);
      }
    }
      
    /* --- Here I am doing Ma*stk + Mb * Su + Mc * S0 + 
           (gam * dS0 - theta * dSu) * dtau / 3.0 to compute the 
           right-hand term
           --                                      ------------------ */
    
    memset(V0, 0, 4*sizeof(double));
    
    for(i = 0;  i < 4;  i++){
      for(j = 0;  j < 4;  j++){
	V0[i] += Ma[i][j] * I[j][k-dk] + Mb[i][j] * Su[j] +
	  Mc[i][j] * S0[j];
      }
      V0[i] += dt03 * (gamma * dS0[i] - theta * dSu[i]);
    }
    /* --- Solve linear system to get the intensity -- -------------- */
      
    SIMD_MatInv(Md[0]);   // Invert Md
    m4v(Md, V0, V1);      // Multiply Md^-1 * V0

    for(i=0;i<4;i++) I[i][k] = V1[i];
      
      
    /* --- Shift values for next depth --          ------------------ */
      
    memcpy(Su,   S0, 4*sizeof(double));
    memcpy(S0,   Sd, 4*sizeof(double));
    memcpy(dSu, dS0, 4*sizeof(double));
      
    memcpy(Ku[0],   K0[0], 16*sizeof(double));
    memcpy(K0[0],   Kd[0], 16*sizeof(double));
    memcpy(dKu[0], dK0[0], 16*sizeof(double));
      
    dtau_uw = dtau_dw;
    dsup    = dsdn;
    dchi_up = dchi_c;
    dchi_c  = dchi_dn;      
  }
      
  /* --- Linear integration in the last interval -- ----------------- */
  
  k = k_end;
  dtau_uw = 0.5*imu * (chi[k] + chi[k-dk]) *
    fabs(geometry.height[k] - geometry.height[k-dk]);
  w3(dtau_uw, w);

  /* --- dSu is defined negative in Han's implementation ------------ */
  
  for (n = 0;  n < 4;  n++)
    V0[n] = w[0]*S[n][k] + w[1] * -dSu[n];
  
  if (Psi) Psi[k] = w[0] - w[1] / dtau_uw;
      
  for (n = 0;  n < 4;  n++) {
    for (m = 0;  m < 4;  m++) {
      A[n][m]  = -w[1]/dtau_uw * Ku[n][m];
      Md[n][m] = (w[0] - w[1]/dtau_uw) * K0[n][m];
    }
    A[n][n]  = 1.0 - w[0];
    Md[n][n] = 1.0;
  }
      
  for (n = 0;  n < 4;  n++) 
    for (m = 0;  m < 4;  m++) 
      V0[n] += A[n][m] * I[m][k-dk];

  /* --- Solve linear system --                    ------------------ */
  
  SIMD_MatInv(Md[0]); // Invert Md
  m4v(Md,V0,V1);      // Multiply Md^-1 * V0
  
  for (n = 0;  n < 4;  n++) I[n][k] = V1[n];
}
/* ------- end ------------------------- Piece_Stokes_Bezier3_1D.c -- */


/* ------- begin ----------------------- Piecewise_Bezier3_1D.c ----- */

void Piecewise_Bezier3_1D(int nspect, int mu, bool_t to_obs,
			  double *chi, double *S, double *I, double *Psi)
{
  
  /* --- Cubic Bezier solver for unpolarized light
         Coded by J. de la Cruz Rodriguez (ISP-SU 2017)

         Reference:
         J. de la Cruz Rodriguez & N. Piskunov (2013), Auer (2003)

         Comments: 
           JdlCR: We only check that the control points of the opacity
	          and source function are always above zero to avoid having
	          a negative interpolant.
         --                                            -------------- */
  
  register int k;
  const char routineName[] = "Piecewise_Bezier3_1D";

  int    k_start, k_end, dk, Ndep = geometry.Ndep;
  double dtau_uw, dtau_dw, dS_uw, I_upw, c1, c2, w[3], zmu, Bnu[2];
  double dsup, dsdn, dt03, eps=0, alpha=0, beta=0, gamma=0, theta=0;
  double dS_up, dS_c, dchi_up, dchi_c, dchi_dn, dsdn2;

  zmu = 1.0 / geometry.muz[mu];

  /* --- Distinguish between rays going from BOTTOM to TOP
         (to_obs == TRUE), and vice versa --           -------------- */

  if (to_obs) {
    dk      = -1;
    k_start = Ndep-1;
    k_end   = 0;
  } else {
    dk      = 1;
    k_start = 0;
    k_end   = Ndep-1;
  }
  
  dtau_uw = 0.5 * zmu * (chi[k_start] + chi[k_start+dk]) *
    fabs(geometry.height[k_start] - geometry.height[k_start+dk]);

  /* --- Boundary conditions --                        -------------- */

  if (to_obs) {
    switch (geometry.vboundary[BOTTOM]) {
    case ZERO:
      I_upw = 0.0;
      break;
    case THERMALIZED:
      Planck(2, &atmos.T[Ndep-2], spectrum.lambda[nspect], Bnu);
      I_upw = Bnu[1] - (Bnu[0] - Bnu[1]) / dtau_uw;
      break;
    case IRRADIATED:
      I_upw = geometry.Ibottom[nspect][mu];
      break;
    case REFLECTIVE:
      sprintf(messageStr, "Boundary condition not implemented: %d",
	      geometry.vboundary[BOTTOM]);
      Error(ERROR_LEVEL_2, routineName, messageStr);
    }
  } else {
    switch (geometry.vboundary[TOP]) {
    case ZERO:
      I_upw = 0.0;
      break;
    case THERMALIZED:
      Planck(2, &atmos.T[0], spectrum.lambda[nspect], Bnu);
      I_upw = Bnu[0] - (Bnu[1] - Bnu[0]) / dtau_uw;
      break;
    case IRRADIATED:
      I_upw = geometry.Itop[nspect][mu];
      break;
    case REFLECTIVE:
      sprintf(messageStr, "Boundary condition not implemented: %d",
	      geometry.vboundary[BOTTOM]);
      Error(ERROR_LEVEL_2, routineName, messageStr);
    }
  }
  
  I[k_start] = I_upw;
  if (Psi) Psi[k_start] = 0.0;
  
  /* --- Set variables for first iteration to allow simple 
         shift for all next iterations --              -------------- */

  k = k_start+dk;
  dsup = fabs(geometry.height[k] - geometry.height[k-dk]) * zmu;
  dsdn = fabs(geometry.height[k+dk] - geometry.height[k]) * zmu;
  dchi_up = (chi[k] - chi[k-dk]) / dsup;
  
  /* --- dchi/ds at central point --                   -------------- */

  dchi_c = cent_deriv(dsup, dsdn, chi[k-dk], chi[k], chi[k+dk]);
  
  /* --- upwind path_length (Bezier3 integration) --   -------------- */

  c1 = MAX(chi[k]    - (dsup/3.0) * dchi_c,   0.0);
  c2 = MAX(chi[k-dk] + (dsup/3.0) * dchi_up,  0.0);
  dtau_uw = dsup * (chi[k] + chi[k-dk] + c1 + c2) * 0.25;

  /* dS/dtau at upwind point */

  dS_up = (S[k] - S[k-dk]) / dtau_uw;

  /* --- Solve transfer along ray --                   -------------- */

  for (k = k_start+dk;  k != k_end+dk;  k += dk) {
    
    if (k != k_end) {

      /* --- Downwind path length --                   -------------- */
      
       dsdn = fabs(geometry.height[k+dk] - geometry.height[k]) * zmu;
       
       /* --- dchi/ds at downwind point --             -------------- */
       
       if (fabs(k - k_end) > 1) {
	 dsdn2=fabs(geometry.height[k+2*dk] -
		    geometry.height[k+dk]) * zmu;
	 dchi_dn = cent_deriv(dsdn,dsdn2,chi[k],chi[k+dk],chi[k+2*dk]);       
       } else {
	 dchi_dn=(chi[k+dk]-chi[k])/dsdn;
       }
       
       /* --- Make sure that c1 and c2 don't go below zero -- ------- */
    
       c1 = MAX(chi[k]    + (dsdn/3.0) * dchi_c,  0.0);
       c2 = MAX(chi[k+dk] - (dsdn/3.0) * dchi_dn, 0.0);

       /* --- Downwind optical path length --          -------------- */

       dtau_dw = dsdn * (chi[k] + chi[k+dk] + c1 + c2) * 0.25;
       dt03    = dtau_uw / 3.0;
      
      /* --- Compute interpolation parameters --       -------------- */
       
       Bezier3_coeffs(dtau_uw, &beta, &alpha, &theta, &gamma, &eps);
       
       /* --- dS/dt at central point --                -------------- */
       
       dS_c = cent_deriv(dtau_uw, dtau_dw, S[k-dk], S[k], S[k+dk]);

       /* --- Source function control points --        -------------- */
       
       c1 = MAX(S[k]    - dt03 * dS_c , 0.0);
       c2 = MAX(S[k-dk] + dt03 * dS_up, 0.0);       
     
       /* --- Solve integral in this interval --       -------------- */
       
       I[k]= I_upw*eps + alpha*S[k] + beta*S[k-dk] +
	 gamma * c1 + theta * c2; 

       /* --- Diagonal operator --                     -------------- */

       if (Psi) Psi[k] = alpha + gamma;
       
    } else { 
      
      /* --- Piecewise linear integration at end of ray -- ---------- */
      
      dtau_uw = 0.5 * zmu * (chi[k] + chi[k-dk]) *
	fabs(geometry.height[k] - geometry.height[k-dk]);
      
      /* --- Defined negative in Han's implementation -- ------------ */
      
      dS_uw = -(S[k] - S[k-dk]) / dtau_uw;
      w3(dtau_uw, w);
      
      I[k] = (1.0 - w[0])*I_upw + w[0]*S[k] + w[1]*dS_uw;

      /* --- Diagonal operator --                      -------------- */
      
      if (Psi) Psi[k] = w[0] - w[1] / dtau_uw;
    }
    
    /* --- Re-use downwind quantities for next upwind position -- --- */
    
    I_upw = I[k];
    dsup=dsdn;
    dchi_up=dchi_c;
    dchi_c=dchi_dn;
    dtau_uw=dtau_dw;
    dS_up = dS_c;
  }
}
/* ------- end ---------------------------- Piecewise_Bezier3_1D.c -- */
