/* ------- file: -------------------------- bezier_aux.c ------------

   Auxiliary routines for Cubic DELO-Bezier (polarized) and
   cubic short-char Bezier solvers.
   
   References: de la Cruz Rodriguez & Piskunov (2013), Auer (2003)
               (Derivatives) Fritsch & Butland (1984),
	       
   Coded by J. de la Cruz Rodriguez (ISP-SU 2017)

   Modifications:
           2017-03-12, JdlCR: Created!

   Adapted for RH:
           Han Uitenbroek

       Last modified: Tue Apr 28 14:30:54 2020 --

       --------------------------                      ----------RH-- */


#include <math.h>
#include <string.h>    // memcpy, memset

#include "bezier.h"


/* --- Global variables --                             -------------- */


/* ------- begin -------------------------- solveLiearFast.c -------- */

void solveLinearFast(double A[4][4], double B[4])
{
  
  // --- Simple Gaussian elimination with partial pivoting --- //
  
  // A is the matrix with coeffs that multiply X. B is the right hand
  // term (on input). The result is overwritten into B. All operations
  // are done in-place
  // It uses the macro "swap" to swap the values of two rows when a pivot
  // is found

  // The algorithm is a pretty standard Gauss-Jordan algorithm, here optimized a
  // bit for performance.
  
  register int i,j,k;
  int maxrow,swapme;
  double maxel,tmp;
  
  
  for (i=0; i<4; i++) {

    // Find pivot

    maxel = fabs(A[i][i]);
    maxrow = i, swapme = 0;
    
    for (k=i+1; k<4; k++){
      tmp = fabs(A[k][i]);
      if(tmp > maxel){
	maxel = tmp;
	maxrow = k, swapme = 1;
      }
    }
    
    // swap
    if(swapme){
      for (k=i; k<4;k++) SWAP(A[maxrow][k],A[i][k],tmp);  
      SWAP(B[maxrow],B[i],tmp); 
    }

    // Set to zero relevant columns
    for (k=i+1; k<4; k++){
      tmp = -A[k][i]/A[i][i];
      A[k][i] = 0.0;
      for ( j=i+1; j<4; j++) {
	A[k][j] += tmp * A[i][j];
      }
      B[k] += tmp*B[i];
    }
  }

  
  // Solve upper triagonal system and store in B
  
  for (i=3; i>=0; i--) {
    B[i] /= A[i][i];
    for ( k=i-1;k>=0; k--) {
      B[k] -= A[k][i] * B[i];
    }
  }
  
}

/* ------- end ---------------------------- solveLinearFast.c ------- */


/* ------- begin -------------------------- cent_deriv.c ------------ */

inline double cent_deriv(double odx,double dx,
			 double yu,double y0, double yd)
{
  /* --- Derivatives from Steffen (1990) --- */
  
  const double S0 = (yd - y0) / dx;
  const double Su = (y0 - yu) / odx;
  const double P0 = fabs((Su*dx + S0*odx) / (odx+dx)) * 0.5;
  
  return (signFortran(S0) + signFortran(Su)) *
    fmin(fabs(Su), fmin(fabs(S0), P0));
}
/* ------- end ---------------------------- cent_deriv.c ------------ */


/* ------- begin -------------------------- cent_deriv_mat.c -------- */

inline void cent_deriv_mat(double wprime[4][4], double dsup, double dsdn,
			   double chiup[4][4], double chic[4][4],
			   double chidn[4][4])
{
  register int i,j;
  
  for(j = 0;  j<4;  j++)
    for(i = 0;  i < 4;  i++)
      wprime[j][i] = cent_deriv(dsup, dsdn, chiup[j][i], chic[j][i],
				chidn[j][i]);
}
/* ------- end ---------------------------- cent_deriv_mat.c -------- */


/* ------- begin -------------------------- cent_deriv_vec.c -------- */

inline void cent_deriv_vec(double wprime[4], double dsup, double dsdn,
		    double chiup[4], double chic[4], double chidn[4])
{
  register int i;
  
  for(i = 0;  i < 4;  i++)
    wprime[i] = cent_deriv(dsup, dsdn, chiup[i], chic[i], chidn[i]);
  
}
/* ------- end ---------------------------- cent_deriv_vec.c -------- */


/* ------- begin -------------------------- m4m.c ------------------- */


inline void m4m(double a[4][4], double b[4][4], double c[4][4])
{

  /* --- Matrix multiplication --                      -------------- */

  register int i, j, k;
  memset(&c[0][0],0,sizeof(double)*16);
  for(j = 0; j<4; j++)
    for(i = 0; i<4; i++)
      for(k = 0; k<4; k++)
	c[j][i] += a[k][i]*b[j][k]; 
}
/* ------- end ---------------------------- m4m.c ------------------- */


/* ------- begin -------------------------- m4v.c ------------------- */

/* --- Matrix/vector multiplication.
       We use matrix as float to be able to use Intel's
       matrix inversion --                         ------------------ */

inline void m4v(float a[4][4], double b[4], double c[4])
{
  register int k, i;
  memset(&c[0],0,sizeof(double)*4);
  for(i = 0; i<4; i++)
    for(k = 0; k<4; k++)
      c[i] += ((double)a[i][k]) * b[k];
}
/* ------- end ---------------------------- m4v.c ------------------- */


/* ------- begin -------------------------- Svec.c ------------------ */

inline void Svec(int k, double **S, double *Sf)
{
  /* --- Extracts the Source vector at depth-points k -- ------------ */

  Sf[0] = S[0][k], Sf[1] = S[1][k], Sf[2] = S[2][k], Sf[3] = S[3][k];
}
/* ------- end ---------------------------- Svec.c ------------------ */


/* ------- begin -------------------------- Bezier3_coeffs ---------- */

inline void Bezier3_coeffs(double dt, double *alpha, double *beta,
			   double *gamma, double *theta, double *eps)
{
  /* ---
     Integration coeffs. for cubic Bezier interpolants
     Use Taylor expansion if dtau is small
     Note:
        alpha -> Su
	beta  -> Sc
	Gamma -> Cc
	theta -> Cu
     
     --- */
  
  double dt2 = dt*dt, dt3 = dt2 * dt;//,  dt4;
  if(dt < 0.05){
    *eps =   1.0 - dt + 0.50 * dt2 - dt3 *  0.166666666666666666666666667;
    *alpha = 0.25 * dt - 0.20 * dt2 +
      dt3 * 0.083333333333333333333333333;// - dt4 / 840.0;
    *beta  = 0.25 * dt - 0.05 * dt2 +
      dt3 * 0.008333333333333333333333333;//  - dt4 / 42.0;
    *gamma = 0.25 * dt - 0.15 * dt2 + dt3 * 0.050;// - dt4 / 210.0;
    *theta = 0.25 * dt - 0.10 * dt2 + dt3 * 0.025;// - dt4 / 84.0;
    return;
  }else{
    if((dt > 30.)){
      *eps = 0.0;
      *alpha = 6.0 / dt3;
      *beta = 1.0 - (*alpha)*(1.0-dt) - 3.0/dt;
      *gamma = *alpha*(dt-3.0);
      *theta = 3.0*(*alpha + (dt-4.0)/dt2);
      return;
    }else{
      *eps = exp(-dt);
      *alpha = -(-6.0+(6.0+6.0*dt+3.0*dt2+dt3)*(*eps))/dt3;
      *beta  = (-6.0 + dt*(6.0+dt*(dt-3.0)) +6.0*(*eps)) / dt3;
      *gamma = 3.0 * (2.0*dt-6.0 + *eps*(6.0+dt*(dt+4.0))) / dt3;
      *theta = 3.0 * (-2.0*(*eps)*(dt+3.0) +6.0+(dt-4.0)*dt) / dt3;
      return;
    }
  }
}
/* ------- end ---------------------------- Bezier3_coeffs ---------- */


/* ------- begin -------------------------- m4inv.c ----------------- */

void m4inv(double MI[4][4])
{

  /* --- In-place Shipley-Coleman matrix inversion
         Fast, but ... how accurate??
         Pivoting is always done in the diagonal.
         Copied here just in case the SIMD matrix inversion 
         gives troubles. --                            -------------- */
  
  register int k, i, j;
  
  for (k = 0;  k < 4;  k++){

    /* --- The pivot element --                        -------------- */
    
    MI[k][k] = -1.0 / MI[k][k];

    /* --- The pivot column --                         -------------- */
    
    for(i = 0;  i < 4;  ++i) if(i != k) MI[i][k]*=MI[k][k];

    /* --- Elements not in a pivot row or column --    -------------- */
    
    for(i = 0;  i < 4;  ++i) {
      if(i != k)
	for(j = 0;  j < 4;  ++j)
	  if(j != k)
	    MI[i][j] += MI[i][k] * MI[k][j];
    }
    /* --- Elements in a pivot row --                  -------------- */
    
    for(i = 0;  i < 4;  ++i) {
      if(i != k)
	MI[k][i] *= MI[k][k];
    }
  }
  
  for(i = 0;  i < 4;  ++i) {
    for(j = 0;  j < 4;  ++j) MI[i][j] = -MI[i][j];
  }
  return;
}
/* ------- end ---------------------------- m4inv.c ----------------- */
