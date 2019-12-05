/* ------- file: -------------------------- piecewise_1D.c ----------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Thu May 31 13:26:58 2018 --

       --------------------------                      ----------RH-- */

/* --- Piecewise integration of transfer equation in one dimension.

    -- For boundary condition THERMALIZED use the relation

          I ~= B - \mu dB / d\tau,

       where \tau is taken in the direction of the ray (i.e. NOT in the
       sense of optical depth, e.g. see Mihalas, 1978, p. 51).
       --                                              -------------- */


#include <math.h>
#include <stdlib.h>

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


/* ------- begin -------------------------- Piecewise_Linear_1D.c --- */

void Piecewise_Linear_1D(int nspect, int mu, bool_t to_obs,
			 double *chi, double *S, double *I, double *Psi)
{
  const char routineName[] = "Piecewise_Linear_1D";
  register int k;

  int    k_start, k_end, dk, Ndep = geometry.Ndep;
  double dtau_uw, dtau_dw, dS_uw, dS_dw, I_uw, w[2], zmu, Bnu[2];

  zmu = 0.5 / geometry.muz[mu];

  /* --- Distinguish between rays going from BOTTOM to TOP
         (to_obs == TRUE), and vice versa --      -------------- */

  if (to_obs) {
    dk      = -1;
    k_start = Ndep-1;
    k_end   = 0;
  } else {
    dk      = 1;
    k_start = 0;
    k_end   = Ndep-1;
  }
  dtau_uw = zmu * (chi[k_start] + chi[k_start+dk]) *
    fabs(geometry.height[k_start] - geometry.height[k_start+dk]);
  dS_uw = (S[k_start] - S[k_start+dk]) / dtau_uw;

  /* --- Boundary conditions --                        -------------- */

  if (to_obs) {
    switch (geometry.vboundary[BOTTOM]) {
    case ZERO:
      I_uw = 0.0;
      break;
    case THERMALIZED:
      Planck(2, &atmos.T[Ndep-2], spectrum.lambda[nspect], Bnu);
      I_uw = Bnu[1] - (Bnu[0] - Bnu[1]) / dtau_uw;
      break;
    case IRRADIATED:
      I_uw = geometry.Ibottom[nspect][mu];
      break;
    case REFLECTIVE:
      sprintf(messageStr, "Boundary condition not implemented: %d",
	      geometry.vboundary[BOTTOM]);
      Error(ERROR_LEVEL_2, routineName, messageStr);
    }
  } else {
    switch (geometry.vboundary[TOP]) {
    case ZERO:
      I_uw = 0.0;
      break;
    case THERMALIZED:
      Planck(2, &atmos.T[0], spectrum.lambda[nspect], Bnu);
      I_uw = Bnu[0] - (Bnu[1] - Bnu[0]) / dtau_uw;
      break;
    case IRRADIATED:
      I_uw = geometry.Itop[nspect][mu];
      break;
    case REFLECTIVE:
      sprintf(messageStr, "Boundary condition not implemented: %d",
	      geometry.vboundary[BOTTOM]);
      Error(ERROR_LEVEL_2, routineName, messageStr);
    }
  }

  I[k_start] = I_uw;
  if (Psi) Psi[k_start] = 0.0;

  /* --- Solve transfer along ray --                   -------------- */

  for (k = k_start+dk;  k != k_end;  k += dk) {
    w2(dtau_uw, w);

    I[k] = (1.0 - w[0])*I_uw + w[0]*S[k] + w[1]*dS_uw;
    if (Psi) Psi[k] = w[0] - w[1] / dtau_uw;

    dtau_uw = zmu * (chi[k] + chi[k+dk]) *
      fabs(geometry.height[k] - geometry.height[k+dk]);
    dS_uw   = (S[k] - S[k+dk]) / dtau_uw;

    I_uw = I[k];
  }
  I[k_end] = (1.0 - w[0])*I_uw + w[0]*S[k_end] + w[1]*dS_uw;
  if (Psi) Psi[k_end] = w[0] - w[1] / dtau_uw;
}
/* ------- end ---------------------------- Piecewise_Linear_1D.c --- */


/* ------- begin -------------------------- Piecewise_1D.c ---------- */

void Piecewise_1D(int nspect, int mu, bool_t to_obs,
		  double *chi, double *S, double *I, double *Psi)
{
  const char routineName[] = "Piecewise_1D";
  register int k;

  int    k_start, k_end, dk, Ndep = geometry.Ndep;
  double dtau_uw, dtau_dw, dS_uw, I_uw, dS_dw, c1, c2, w[3],
         zmu, Bnu[2];

  zmu = 0.5 / geometry.muz[mu];

  /* --- Distinguish between rays going from BOTTOM to TOP
         (to_obs == TRUE), and vice versa --      -------------- */

  if (to_obs) {
    dk      = -1;
    k_start = Ndep-1;
    k_end   = 0;
  } else {
    dk      = 1;
    k_start = 0;
    k_end   = Ndep-1;
  }
  
  dtau_uw = zmu * (chi[k_start] + chi[k_start+dk]) *
    fabs(geometry.height[k_start] - geometry.height[k_start+dk]);

  /* --- Boundary conditions --                        -------------- */

  if (to_obs) {
    switch (geometry.vboundary[BOTTOM]) {
    case ZERO:
      I_uw = 0.0;
      break;
    case THERMALIZED:
      Planck(2, &atmos.T[Ndep-2], spectrum.lambda[nspect], Bnu);
      I_uw = Bnu[1] - (Bnu[0] - Bnu[1]) / dtau_uw;
      break;
    case IRRADIATED:
      I_uw = geometry.Ibottom[nspect][mu];
      break;
    case REFLECTIVE:
      sprintf(messageStr, "Boundary condition not implemented: %d",
	      geometry.vboundary[BOTTOM]);
      Error(ERROR_LEVEL_2, routineName, messageStr);
    }
  } else {
    switch (geometry.vboundary[TOP]) {
    case ZERO:
      I_uw = 0.0;
      break;
    case THERMALIZED:
      Planck(2, &atmos.T[0], spectrum.lambda[nspect], Bnu);
      I_uw = Bnu[0] - (Bnu[1] - Bnu[0]) / dtau_uw;
      break;
    case IRRADIATED:
      I_uw = geometry.Itop[nspect][mu];
      break;
    case REFLECTIVE:
      sprintf(messageStr, "Boundary condition not implemented: %d",
	      geometry.vboundary[BOTTOM]);
      Error(ERROR_LEVEL_2, routineName, messageStr);
    }
  }

  I[k_start] = I_uw;
  if (Psi) Psi[k_start] = 0.0;

  /* --- Solve transfer along ray --                   -------------- */

  dS_uw = (S[k_start] - S[k_start+dk]) / dtau_uw;

  for (k = k_start+dk;  k != k_end+dk;  k += dk) {
    w3(dtau_uw, w);

    if (k != k_end) {

      /* --- Piecewise quadratic here --               -------------- */ 

      dtau_dw = zmu * (chi[k] + chi[k+dk]) *
	fabs(geometry.height[k] - geometry.height[k+dk]);
      dS_dw   = (S[k] - S[k+dk]) / dtau_dw;

      c1 = (dS_uw*dtau_dw + dS_dw*dtau_uw);
      c2 = (dS_uw - dS_dw);
      I[k] = (1.0 - w[0])*I_uw + w[0]*S[k] +
	(w[1]*c1 + w[2]*c2) / (dtau_uw + dtau_dw);

      /* --- Try piecewise linear if quadratic gives negative
             monochromatic intensity --                -------------- */ 

      if (I[k] < 0.0) {
	c1   = dS_uw;
	I[k] = (1.0 - w[0])*I_uw + w[0]*S[k] + w[1]*c1;

	if (Psi) Psi[k] = w[0] - w[1]/dtau_uw;
      } else {
	if (Psi) {
	  c1 = dtau_uw - dtau_dw;
	  Psi[k] = w[0] + (w[1]*c1 - w[2]) / (dtau_uw * dtau_dw);
	}
      }
    } else {
	
      /* --- Piecewise linear integration at end of ray -- ---------- */

      I[k] = (1.0 - w[0])*I_uw + w[0]*S[k] + w[1]*dS_uw;
      if (Psi) Psi[k] = w[0] - w[1] / dtau_uw;
    }
    I_uw = I[k];

    /* --- Re-use downwind quantities for next upwind position -- --- */

    dS_uw   = dS_dw;
    dtau_uw = dtau_dw;
  }
}
/* ------- end ---------------------------- Piecewise_1D.c ---------- */
