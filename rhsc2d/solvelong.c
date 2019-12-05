/* ------- file: -------------------------- solvelong.c -------------

       Version:       rh2.0, 2-D Cartesian
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Thu May 31 11:17:51 2018 --

       --------------------------                      ----------RH-- */

/* --- Piecewise integration of source function for 
       long characteristics. --                        -------------- */


#include <math.h>
#include <stdlib.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "spectrum.h"
#include "geometry.h"
#include "inputs.h"


/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

extern InputData input;


/* ------- begin -------------------------- SolveLong.c ------------- */

/* --- Solves transfer equation along a long characteristic specified in
       LongChar structure pointed to by lc. Used when, for periodic
       boundary conditions, no shortcharacteristic crosses a horizontal
       grid line in the upwind direction. --           -------------- */


double SolveLong(LongChar *lc, int local, double *chi, double *S, double *I)
{
  register int ls;

  bool_t monotonic = TRUE;
  double c1, c2, dtau_dw, dS_uw, dS_dw, w[3], chi_loc, S_loc,
         dtau_uw, chi_uw, S_uw, I_uwp, I_uw, chi_dw, S_dw,
         U[3], psi_uw, psi_dw, Sc, psi_0;

  /* --- The first point is the intersection with the horizontal
         grid line --                                  -------------- */

  S_uw   = Quadr(S,   &lc->stencil[0], UPWIND, monotonic); 
  chi_uw = Quadr(chi, &lc->stencil[0], UPWIND, monotonic);
  I_uwp  = Quadr(I,   &lc->stencil[0], UPWIND, monotonic);

  S_loc   = Quadr(S,   &lc->stencil[1], UPWIND, monotonic);
  chi_loc = Quadr(chi, &lc->stencil[1], UPWIND, monotonic);
  dtau_uw = 0.5 * (chi_uw + chi_loc) * lc->stencil[0].ds[UPWIND];

  /* --- When lc->Nst == 2 we need I_uw = I_uwp --     -------------- */

  I_uw = I_uwp;

  /* --- Go through the long characteristic section by section -- --- */

  for (ls = 2;  ls < lc->Nlc;  ls++) {
    if (ls == lc->Nlc-1) {

      /* --- The last point is the the endpoint for which the non-local
             contribution is needed. --                -------------- */

      S_dw   = S[local];
      chi_dw = chi[local];
    } else {
      S_dw   = Quadr(S,   &lc->stencil[ls], UPWIND, monotonic);
      chi_dw = Quadr(chi, &lc->stencil[ls], UPWIND, monotonic);
    }
    dtau_dw = 0.5 * (chi_loc + chi_dw) * lc->stencil[ls-1].ds[UPWIND];

    /* --- Parabolic interpolation of source function -- ------------ */

    w3(dtau_uw, w);

    dS_uw = (S_uw - S_loc) / dtau_uw;
    dS_dw = (S_loc - S_dw) / dtau_dw;

    c1 = dS_uw*dtau_dw + dS_dw*dtau_uw;
    c2 = dS_uw - dS_dw;
    I_uw = I_uwp*(1.0 - w[0]) + w[0]*S_loc +
      (w[1]*c1 + w[2]*c2) / (dtau_uw + dtau_dw);

    /* --- Try piecewise linear if quadratic gives negative
       monochromatic intensity --                -------------- */ 

    if (I_uw < 0.0) {
      c1   = dS_uw;
      I_uw = (1.0 - w[0])*I_uwp + w[0]*S_loc + w[1]*c1;
    }

    I_uwp = I_uw;
    chi_uw  = chi_loc;  S_uw  = S_loc;
    chi_loc = chi_dw;   S_loc = S_dw;
    dtau_uw = dtau_dw;
  }

  return I_uw;
}
/* ------- end ---------------------------- SolveLong.c ------------- */

/* ------- begin -------------------------- SolveLongStokes.c ------- */

void SolveLongStokes(int nspect, LongChar *lc, int local, double *chi,
		     double **S, double **I, double *I_uw)
{
  register int ls, n, m;

  bool_t monotonic = TRUE;
  double c1, c2, dtau_dw, dS_uw[4], dS_dw[4], w[3],
         chi_loc, S_loc[4], dtau_uw, chi_uw, S_uw[4], chi_dw, S_dw[4],
         P[4], Q[4][4], **R, K[4][4], K_uw[4][4];

  R = matrix_double(4, 4);

  /* --- The first point is the intersection with the horizontal
         grid line --                                  -------------- */

  chi_uw = Quadr(chi, &lc->stencil[0], UPWIND, monotonic);
  StokesK_2D(nspect, &lc->stencil[0], UPWIND, monotonic, chi_uw, K_uw);

  for (n = 0;  n < 4;  n++) {
    S_uw[n]  = Quadr(S[n], &lc->stencil[0], UPWIND, monotonic); 
    I_uw[n]  = Quadr(I[n], &lc->stencil[0], UPWIND, monotonic);
    S_loc[n] = Quadr(S[n], &lc->stencil[1], UPWIND, monotonic);
  }
  chi_loc = Quadr(chi, &lc->stencil[1], UPWIND, monotonic);
  StokesK_2D(nspect, &lc->stencil[1], UPWIND, monotonic, chi_loc, K);

  dtau_uw = 0.5 * (chi_uw + chi_loc) * lc->stencil[0].ds[UPWIND];

  for (ls = 2;  ls < lc->Nlc;  ls++) {
    if (ls == lc->Nlc-1) {

      /* --- The last point is the the endpoint for which the non-local
             contribution is needed. --                -------------- */

      chi_dw = chi[local];
      for (n = 0;  n < 4;  n++) S_dw[n] = S[n][local];
    } else {
      chi_dw = Quadr(chi, &lc->stencil[ls], UPWIND, monotonic);
      for (n = 0;  n < 4;  n++)
	S_dw[n] = Quadr(S[n], &lc->stencil[ls], UPWIND, monotonic);
    }
    dtau_dw = 0.5 * (chi_loc + chi_dw) * lc->stencil[ls-1].ds[UPWIND];
    w3(dtau_uw, w);

    for (n = 0;  n < 4;  n++) {
      dS_uw[n] = (S_uw[n] - S_loc[n]) / dtau_uw;
      dS_dw[n] = (S_loc[n] - S_dw[n]) / dtau_dw;
      c1 = dS_uw[n]*dtau_dw + dS_dw[n]*dtau_uw;
      c2 = dS_uw[n] - dS_dw[n];
      P[n] = w[0]*S_loc[n] + (w[1]*c1 + w[2]*c2) / (dtau_uw + dtau_dw);
    }
    for (n = 0;  n < 4;  n++) {
      for (m = 0;  m < 4;  m++) {
	Q[n][m] = -w[1]/dtau_uw * K_uw[n][m];
	R[n][m] = (w[0] - w[1]/dtau_uw) * K[n][m];
      }
      Q[n][n] = 1.0 - w[0];
      R[n][n] = 1.0;
    }
    for (n = 0;  n < 4;  n++) {
      for (m = 0;  m < 4;  m++) 
	P[n] += Q[n][m] * I_uw[m];
    }
    /* --- Solve linear equations for I --             -------------- */

    SolveLinearEq(4, R, P, TRUE);

    /* --- Store results for the upwind Stokes vector -- ------------ */
    
    for (n = 0;  n < 4;  n++) I_uw[n] = P[n];

    /* --- Reuse upwind quantities --                  -------------- */

    if (ls < lc->Nlc-1) {
      chi_uw  = chi_loc;
      chi_loc = chi_dw;
      dtau_uw = dtau_dw;

      for (n = 0;  n < 4;  n++) {
	S_uw[n]  = S_loc[n];
	S_loc[n] = S_dw[n];
	for (m = 0;  m < 4;  m++) K_uw[n][m] = K[n][m];
      }
      StokesK_2D(nspect, &lc->stencil[ls],
		 UPWIND, monotonic, chi_dw, K);
    }
  }
  freeMatrix((void **) R);
}
/* ------- end ---------------------------- SolveLongStokes.c ------- */
