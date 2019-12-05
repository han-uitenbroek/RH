/* ------- file: -------------------------- feautrier.c -------------

       Version:       rh2.0, 1-D spherically symmetric
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Wed May 23 10:52:47 2018 --

       --------------------------                      ----------RH-- */

/* --- Evaluate sum of monochromatic intensities (I^+ + I^-)/2 and
       return emergent intensity I_surface along a ray with given
       optical depth scale tau and source function S.

    -- When FeautrierOrder == FEAUTRIER_HERMITE the fourth order Feautrier method
       suggested by L. Auer (1976, JQSRT 16, 931-938) is used.

    -- Numerical scheme formulated as in Appendix A of:
       G.B. Rybicki & D.G Hummer, A&A 245, 171-181.

       Boundary -- upper: Iminus = r0*Iplus  + h0
                   lower: Iplus  = rN*Iminus + hN
       --                                              -------------- */
 
#include <stdlib.h>
#include <math.h>

#include "rh.h"
#include "error.h"
#include "atom.h"
#include "atmos.h"
#include "geometry.h"
#include "spectrum.h"
#include "error.h"

#define  A_SIXTH  0.166666666667


/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

extern Geometry geometry;
extern Atmosphere atmos;
extern Spectrum spectrum;
extern char messageStr[];


/* ------- begin -------------------------- Feautrier.c ------------- */

double Feautrier(int nspect, int mu, double *chi, double *S,
		 enum FeautrierOrder F_order, double *P, double *Psi)
{
  const char routineName[] = "Feautrier";
  register int k;

  int    Ns;
  double r0, h0, rN, hN, f0, fN, Ak, Ck, tau0 = 0.0, Bnu[2], dtau_mid,
         Iplus, *dtau, *abc, *A1, *C1, *F, *G, *Q, *Stmp, *ztmp;
  Ray   *ray;

  ray = &geometry.rays[mu];
  Ns  = ray->Ns;

  dtau = (double *) malloc(Ns * sizeof(double));
  abc  = (double *) malloc(Ns * sizeof(double));
  Q    = (double *) malloc(Ns * sizeof(double));
  A1   = (double *) malloc(Ns * sizeof(double));
  C1   = (double *) malloc(Ns * sizeof(double));
  F    = (double *) malloc(Ns * sizeof(double));
  G    = (double *) malloc(Ns * sizeof(double));
  ztmp = (double *) malloc(Ns * sizeof(double));
  Stmp = (double *) malloc(Ns * sizeof(double));

  for (k = 0;  k < Ns-1;  k++)
    dtau[k] = 0.5 * (chi[k] + chi[k+1]) * (ray->s[k] - ray->s[k+1]);

  /* --- Upper boundary condition:  I^- = r0*I^+ + h0 --   ---------- */

  switch (geometry.rboundary[TOP]) {
  case ZERO:
    r0 = h0 = 0.0;
    break;
  case IRRADIATED:
    r0 = 0.0;
    h0 = geometry.Itop[nspect];
    break;
    case THERMALIZED:
      sprintf(messageStr, "Boundary condition not implemented: %d",
	      geometry.rboundary[TOP]);
      Error(ERROR_LEVEL_2, routineName, messageStr);
  }
  f0      = (1.0 - r0) / (1.0 + r0);
  abc[0]  = 1.0 + 2.0*f0 / dtau[0];
  C1[0]   = 2.0 / SQ(dtau[0]);
  Stmp[0] = S[0] + 2.0*h0 / ((1.0 + r0)*dtau[0]);
  if (F_order == FEAUTRIER_HERMITE) {
    C1[0]   -= 2.0*A_SIXTH;
    Stmp[0] += 2.0*A_SIXTH * (S[1] - S[0]);
  }
  /* --- Lower boundary condition:  I^+ = rN*I^- + hN.
         If ray is a SHELL_RAY use reflective boundary, otherwise use
         THERMALIZED. --                               -------------- */

  if (ray->type == SHELL_RAY) {
    rN = 1.0;
    hN = 0.0;
  } else {
    rN = 0.0;
    Planck(2, &atmos.T[geometry.Nradius-2], spectrum.lambda[nspect], Bnu);
    hN = Bnu[1] - (Bnu[0] - Bnu[1]) / dtau[geometry.Nradius-2];
  }
  fN = (1.0 - rN) / (1.0 + rN);
  abc[Ns-1]  = 1.0 + 2.0*fN / dtau[Ns-2];
  A1[Ns-1]   = 2.0 / SQ(dtau[Ns-2]);
  Stmp[Ns-1] = S[Ns-1] + 2.0*hN / ((1.0 + rN)*dtau[Ns-2]);
  if (F_order == FEAUTRIER_HERMITE) {
    A1[Ns-1]   -= 2.0*A_SIXTH;
    Stmp[Ns-1] += 2.0*A_SIXTH * (S[Ns-2] - S[Ns-1]);
  }

  for (k = 1;  k < Ns-1;  k++) {
    dtau_mid = 0.5*(dtau[k] + dtau[k-1]);
    A1[k]   = 1.0 / (dtau_mid * dtau[k-1]);
    C1[k]   = 1.0 / (dtau_mid * dtau[k]);
    abc[k]  = 1.0;
    Stmp[k] = S[k];
  }
  if (F_order == FEAUTRIER_HERMITE) {
    for (k = 1;  k < Ns-1;  k++) {
      Ak     = A_SIXTH * (1.0 - 0.5*SQ(dtau[k])*A1[k]);
      Ck     = A_SIXTH * (1.0 - 0.5*SQ(dtau[k-1])*C1[k]);
      A1[k] -= Ak;
      C1[k] -= Ck;
      Stmp[k] += Ak*(S[k-1] - S[k]) + Ck*(S[k+1] - S[k]);
    }
  }
  /* --- Start the elimination --                          ---------- */

  F[0]    = abc[0] / C1[0];
  ztmp[0] = Stmp[0] / (abc[0] + C1[0]);
  for (k = 1;  k < Ns-1;  k++) {
    F[k]    = (abc[k] + A1[k]*F[k-1]/(1.0 + F[k-1])) / C1[k];
    ztmp[k] = (Stmp[k] + A1[k]*ztmp[k-1]) / (C1[k] * (1.0 + F[k]));
  }
  /* --- Now backsubstitution --                          ----------- */

  P[Ns-1] = (Stmp[Ns-1]+ A1[Ns-1]*ztmp[Ns-2]) /
    (abc[Ns-1] + A1[Ns-1]*(F[Ns-2] / (1.0 + F[Ns-2])));
  for (k = Ns-2;  k >= 0;  k--)
    P[k] = P[k+1] / (1.0 + F[k]) + ztmp[k];

  /* --- If necessary evaluate the diagonal operator -- ------------- */

  if (Psi) {
    if (F_order == FEAUTRIER_HERMITE) {
      sprintf(messageStr,
	      "Higher order for diagonal operator not yet implemented");
      Error(ERROR_LEVEL_1, routineName, messageStr);
    }

    G[Ns-1] = abc[Ns-1] / A1[Ns-1];
    for (k = Ns-2;  k >= 1;  k--)
      G[k] = (abc[k] + C1[k]*G[k+1]/(1.0 + G[k+1])) / A1[k];

    Psi[0] = 1.0 / (abc[0] + C1[0]*G[1]/(1.0 + G[1]));
    for (k = 1;  k < Ns-1;  k++)
      Psi[k] = 1.0 / (abc[k] + A1[k]*F[k-1]/(1.0 + F[k-1]) +
		      C1[k]*G[k+1]/(1.0 + G[k+1]));
    Psi[Ns-1] = 1.0 /
      (abc[Ns-1] + A1[Ns-1]*F[Ns-2]/(1.0 + F[Ns-2]));
  }

  free(dtau);  free(abc);  free(Q);   free(A1);
  free(C1);    free(F);    free(G);   free(ztmp);
  free(Stmp);

  /* --- Return emergent intensity --                   ------------- */

  Iplus = (1.0 + f0)*P[0] - h0/(1.0 + r0);

  if (tau0) 
    return (Iplus - S[0])*exp(-tau0) + S[0];
  else
    return Iplus;
}
/* ------- end ---------------------------- Feautrier.c ------------- */
