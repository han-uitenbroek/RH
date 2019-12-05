/* ------- file: -------------------------- piecestokes_1D.c --------

       Version:       rh2.0, 1-D plane-parallel
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Thu May 24 14:21:55 2018 --

       --------------------------                      ----------RH-- */

/* --- Piecewise integration of the coupled Stokes transfer equations.
       Method is quasi-parabolic DELO method.

  See: - D. E. Rees, G. A. Murphy and C. J. Durrant 1989, ApJ 339,
         1093-1106.

       - H. Socas Navarro, J. Trujillo Bueno and B. Ruiz Cobo 2000,
         "Non-LTE Inversion of Stokes Profiles", ApJ 530, 977.

    -- For boundary condition THERMALIZED use the relation

          (I, Q, U, V) ~= (B - \mu dB / d\tau, 0, 0, 0)

       where dB/\tau is taken in the direction of the ray (i.e. NOT in
       the sense of optical depth, e.g. see Mihalas, 1978, p. 51).
       --                                              -------------- */

#include <math.h>

#include "rh.h"
#include "error.h"
#include "atom.h"
#include "atmos.h"
#include "geometry.h"
#include "spectrum.h"


/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

extern Geometry geometry;
extern Atmosphere atmos;
extern Spectrum spectrum;
extern char messageStr[];


/* ------- begin -------------------------- Piece_Stokes_1D.c ------- */

void Piece_Stokes_1D(int nspect, int mu, bool_t to_obs,
		     double *chi_I, double **S, double **I, double *Psi)
{
  const char routineName[] = "Piece_Stokes_1D";
  register int k, n, m;

  int    Ndep = geometry.Ndep, k_start, k_end, dk;
  double dtau_uw, dtau_dw = 0.0, dS_uw[4], dS_dw[4], c1, c2, w[3],
    I_upw[4], zmu, P[4], Bnu[2], Q[4][4], **R, K[4][4], K_upw[4][4];

  zmu = 0.5 / geometry.muz[mu];
  R = matrix_double(4, 4);

  if (to_obs) {
    dk      = -1;
    k_start = Ndep-1;
    k_end   = 0;
  } else {
    dk      = 1;
    k_start = 0;
    k_end   = Ndep-1;
  }
  dtau_uw = zmu * (chi_I[k_start] + chi_I[k_start+dk]) *
    fabs(geometry.height[k_start] - geometry.height[k_start+dk]);

  StokesK(nspect, k_start, chi_I[k_start], K_upw);

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
  for (n = 0;  n < 4;  n++)
    dS_uw[n] = (S[n][k_start] - S[n][k_start+dk]) / dtau_uw;

  for (n = 0;  n < 4;  n++) I[n][k_start] = I_upw[n];
  if (Psi) Psi[k_start] = 0.0;

  for (k = k_start+dk;  k != k_end+dk;  k += dk) {
    w3(dtau_uw, w);
    StokesK(nspect, k, chi_I[k], K);

    if (k != k_end) {
      dtau_dw = zmu * (chi_I[k] + chi_I[k+dk]) *
	fabs(geometry.height[k] - geometry.height[k+dk]);

      for (n = 0;  n < 4;  n++) {
	dS_dw[n] = (S[n][k] - S[n][k+dk]) / dtau_dw;
	c1 = dS_uw[n]*dtau_dw + dS_dw[n]*dtau_uw;
	c2 = dS_uw[n] - dS_dw[n];
	P[n] = w[0]*S[n][k] + (w[1]*c1 + w[2]*c2) / (dtau_uw + dtau_dw);
      }
      if (Psi) {
	c1 = dtau_uw - dtau_dw;
	Psi[k] = w[0] + (w[1]*c1 - w[2]) / (dtau_uw * dtau_dw);
      }
    } else {

      /* --- Piecewise linear integration at end of ray -- ---------- */

      for (n = 0;  n < 4;  n++)	P[n] = w[0]*S[n][k] + w[1]*dS_uw[n];
      if (Psi) Psi[k] = w[0] - w[1] / dtau_uw;
    }

    for (n = 0;  n < 4;  n++) {
      for (m = 0;  m < 4;  m++) {
	Q[n][m] = -w[1]/dtau_uw * K_upw[n][m];
	R[n][m] = (w[0] - w[1]/dtau_uw) * K[n][m];
      }
      Q[n][n] = 1.0 - w[0];
      R[n][n] = 1.0;
    }
    for (n = 0;  n < 4;  n++) {
      for (m = 0;  m < 4;  m++) 
	P[n] += Q[n][m] * I_upw[m];
    }
    /* --- Solve linear equations for I --             -------------- */

    SolveLinearEq(4, R, P, TRUE);

    /* --- Store results for Stokes vector --          -------------- */

    for (n = 0;  n < 4;  n++) I[n][k] = P[n];

    /* --- Re-use downwind quantities for next upwind position -- --- */

    dtau_uw = dtau_dw;
    for (n = 0;  n < 4;  n++) {
      I_upw[n] = I[n][k];
      dS_uw[n] = dS_dw[n];
      for (m = 0;  m < 4;  m++) K_upw[n][m] = K[n][m];
    }
  }
  freeMatrix((void **) R);
}
/* ------- end ---------------------------- Piece_Stokes_1D.c ------- */
