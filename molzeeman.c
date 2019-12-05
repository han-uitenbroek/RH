/* ------- file: -------------------------- molzeeman.c -------------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Thu Jun 30 15:35:10 2011 --

       --------------------------                      ----------RH-- */

/* --- Routines to calculate Zeeman splitting patterns and relative
       strengths of Zeeman components in molecular transitions.
       --                                              -------------- */

#include <stdlib.h>
#include <math.h>

#include "rh.h"
#include "atom.h"
#include "spectrum.h"
#include "error.h"


/* --- Function prototypes --                          -------------- */

double MolZeemanStr(double Ju, double Mu, double Jl, double Ml);
double MolLande_a(double Lambda, double Sigma, double Omega, double J);
double MolLande_b(double Lambda, double S, double N, double J);


/* --- Global variables --                             -------------- */

extern char messageStr[];


/* ------- begin -------------------------- MolZeemanStr.c ---------- */

double MolZeemanStr(double Ju, double Mu, double Jl, double Ml)
{
  const char routineName[] = "MolZeemanStr";

  int    q, dJ;
  double s;

  /* --- Return the strength of Zeeman component (Ju, Mu) -> (Jl, Ml),
         where J and M are the total angular momentum and magnetic
         quantum numbers of the upper and lower level of a Zeeman split
         molecular transition.

    See: - G. Herzberg 1950, in "Spectra of Diatomic Molecules", p.301-3
 
         - S.V. Berdyugina and S.K. Solanki 2002, A&A 385, 701-715

         - A. Schadee 1978, JQSRT 19, 517-531

   Note: Berdyugina & Solanki adhere to the definition that dJ = Jl - Ju

         --                                            -------------- */

  q  = (int) (Ml - Mu);
  dJ = (int) (Jl - Ju);

  switch (dJ) {
  case 1:
    switch (q) {
    case  1: s = (Ju + 1.0 + Mu) * (Ju + 2.0 + Mu) /
	       (2*(Ju + 1.0) * (2*Ju + 1.0) * (2*Ju + 3.0));  break;
    case  0: s = (Ju + 1.0 + Mu) * (Ju + 1.0 - Mu) /
	       ((Ju + 1.0) * (2*Ju + 1.0) * (2*Ju + 3.0));    break;
    case -1: s = (Ju + 1.0 - Mu) * (Ju + 2.0 - Mu) /
	       (2*(Ju + 1.0) * (2*Ju + 1.0) * (2*Ju + 3.0));  break;
    }
    break;

  case 0:
    switch (q) {
    case  1: s = (Ju - Mu) * (Ju + 1.0 + Mu) /
	       (2*Ju * (Ju + 1.0) * (2*Ju + 1.0));  break;
    case  0: s = Mu*Mu /
	       (Ju * (Ju + 1.0) * (2*Ju + 1.0));    break;
    case -1: s = (Ju + Mu) * (Ju + 1.0 - Mu) /
	       (2*Ju * (Ju + 1.0) * (2*Ju + 1.0));  break;
    }
    break;

  case -1:
    switch (q) {
    case  1: s = (Ju - Mu) * (Ju - 1.0 - Mu) /
	       (2*Ju * (2*Ju - 1.0) * (2*Ju + 1.0));  break;
    case  0: s = (Ju + Mu) * (Ju - Mu) /
	       (Ju * (2*Ju - 1.0) * (2*Ju + 1.0));    break;
    case -1: s = (Ju + Mu) * (Ju - 1.0 + Mu) /
	       (2*Ju * (2*Ju - 1.0) * (2*Ju + 1.0));  break;
    }
    break;

  default:
    sprintf(messageStr, "Invalid dJ: %d", dJ);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }   
  /* --- The strengths are normalized already. --      -------------- */

  return s;
}
/* ------- end ---------------------------- MolZeemanStr_b.c -------- */

/* ------- begin -------------------------- MolLande_a.c ------------ */

double MolLande_a(double Lambda, double Sigma, double Omega, double J)
{
  /* --- Lande factor for level in Hund's case a

    See: S.V. Berdyugina and S.K. Solanki 2002, A&A 385, 701-715, Eq. 2
     --                                                -------------- */

  return (Lambda + 2.0*Sigma) * Omega / (J * (J + 1.0));
}
/* ------- end ---------------------------- MolLande_a.c ------------ */

/* ------- begin -------------------------- MolLande_b.c ------------ */

double MolLande_b(double Lambda, double S, double N, double J)
{
  /* --- Lande factor for level in Hund's case b

    See: S.V. Berdyugina and S.K. Solanki 2002, A&A 385, 701-715, Eq. 10
     --                                                -------------- */

  if (Lambda == 0) {
    return 1.0 / (J*(J + 1.0)) *
      (J*(J + 1.0) - N*(N + 1.0) + S*(S + 1.0));
  } else {
    return 1.0 / (J*(J + 1.0)) *
      (Lambda*Lambda / (2*N*(N + 1.0)) *
       (J*(J + 1.0) + N*(N + 1.0) - S*(S + 1.0)) +
       J*(J + 1.0) - N*(N + 1.0) + S*(S + 1.0));
  }
}
/* ------- end ---------------------------- MolLande_b.c ------------ */

/* ------- begin -------------------------- MolLande_eff.c ---------- */

double MolLande_eff(MolecularLine *mrt)
{
  const char routineName[] = "MolLande_eff";
  register double Ml, Mu;

  double g_eff, Ju, Jl, norm, Nl, Nu, shift, strength, gLu, gLl;

  norm  = 0.0;
  g_eff = 0.0;

  Jl = (mrt->gi - 1.0) / 2.0;
  Ju = (mrt->gj - 1.0) / 2.0;

  switch (mrt->Hundi) {
  case CASE_A:
    gLl = MolLande_a(mrt->Lambdai, mrt->Si, mrt->Omegai, Jl);
    break;
  case CASE_B:
    Nl  = Jl - mrt->Si + (mrt->subi - 1);
    gLl = MolLande_b(mrt->Lambdai, mrt->Si, Nl, Jl);
    break;
  default:
    sprintf(messageStr, "Unsupported Hund's case: %d", mrt->Hundi);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }
  switch (mrt->Hundj) {
  case CASE_A:
    gLu = MolLande_a(mrt->Lambdaj, mrt->Sj, mrt->Omegaj, Ju);
    break;
  case CASE_B:
    Nu  = Ju - mrt->Sj + (mrt->subj - 1);
    gLu = MolLande_b(mrt->Lambdaj, mrt->Sj, Nu, Ju);
    break;
  default: 
    sprintf(messageStr, "Unsupported Hund's case: %d", mrt->Hundj);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }

  for (Ml = -Jl;  Ml <= Jl;  Ml++) {
    for (Mu = -Ju;  Mu <= Ju;  Mu++) {
      if ((Ml - Mu) == 1.0) {
	shift    = gLl*Ml - gLu*Mu;
        strength = MolZeemanStr(Ju, Mu, Jl, Ml);

	g_eff += shift * strength;
        norm  += strength;
      }
    }
  }

  return g_eff / norm;
}

/* ------- end ---------------------------- MolLande_eff.c ---------- */

/* ------- begin -------------------------- MolZeeman.c ------------- */

ZeemanMultiplet* MolZeeman(MolecularLine *mrt)
{
  const char routineName[] = "MolZeeman";

  register int n;

  double Jl, Ju, Mu, Ml, norm[3], g_eff, Nl, Nu, gLl, gLu, lambda_air;
  ZeemanMultiplet *zm;

  /* --- Return a pointer to a ZeemanMultiplet structure with all the
         components of a Zeeman split line. The strengths in the line
         are normalized to unity for each of the three possible values
         of q = [-1, 0, 1].

         Convention:
 
          -- q = +1 corresponds to a redshifted \sigma profile
	     (zm->shift > 0). This redshifted profile has
             right-handed circular polarization when the
             magnetic field parallel to the line of sight and
             points towards the observer.

          -- q = 0 corresponds to an unpolarized \pi profile

	 --                                            -------------- */

  zm = (ZeemanMultiplet *) malloc(sizeof(ZeemanMultiplet));

  if (mrt->g_Lande_eff != 0.0) {

    /* --- In case an effective Landee factor has been specified, or
           the the inputs.use_effective_Lande parameter has been set
           in the input --                             -------------- */

    zm->Ncomponent = 3;
    zm->q        = (int *) malloc(3 * sizeof(int));
    zm->strength = (double *) malloc(3 * sizeof(double));
    zm->shift    = (double *) malloc(3 * sizeof(double));

    /* --- Normal Zeeman triplet --                    -------------- */

    for (n = 0;  n < 3;  n++) {
      zm->q[n] = -1 + n;
      zm->strength[n] = 1.0;
      zm->shift[n] = zm->q[n] * mrt->g_Lande_eff;
    }
  } else {

    /* --- Anomalous Zeeman splitting. First, count the number of
           components --                               -------------- */

    Jl = (mrt->gi - 1.0) / 2.0;
    Ju = (mrt->gj - 1.0) / 2.0;

    switch (mrt->Hundi) {
    case CASE_A:
      gLl = MolLande_a(mrt->Lambdai, mrt->Si, mrt->Omegai, Jl);
      break;
    case CASE_B:
      Nl  = Jl - mrt->Si + (mrt->subi - 1);
      gLl = MolLande_b(mrt->Lambdai, mrt->Si, Nl, Jl);
      break;
    default:
      sprintf(messageStr, "Unsupported Hund's case: %d", mrt->Hundi);
      Error(ERROR_LEVEL_2, routineName, messageStr);
    }
    switch (mrt->Hundj) {
    case CASE_A:
      gLu = MolLande_a(mrt->Lambdaj, mrt->Sj, mrt->Omegaj, Ju);
      break;
    case CASE_B:
      Nu  = Ju - mrt->Sj + (mrt->subj - 1);
      gLu = MolLande_b(mrt->Lambdaj, mrt->Sj, Nu, Ju);
      break;
    default: 
      sprintf(messageStr, "Unsupported Hund's case: %d", mrt->Hundj);
      Error(ERROR_LEVEL_2, routineName, messageStr);
    }
    /* --- Determine the number of components --       -------------- */

    zm->Ncomponent = 0;
    for (Ml = -Jl;  Ml <= Jl;  Ml++) {
      for (Mu = -Ju;  Mu <= Ju;  Mu++)
	if (fabs(Mu - Ml) <= 1.0) zm->Ncomponent++;
    }
    zm->q        = (int *) malloc(zm->Ncomponent * sizeof(int));
    zm->strength = (double *) malloc(zm->Ncomponent * sizeof(double));
    zm->shift    = (double *) malloc(zm->Ncomponent * sizeof(double));

    /* --- Fill the structure and normalize the strengths -- -------- */

    for (n = 0;  n < 3;  n++) norm[n] = 0.0;

    n = 0;
    g_eff = 0.0;
    for (Ml = -Jl;  Ml <= Jl;  Ml++) {
      for (Mu = -Ju;  Mu <= Ju;  Mu++) {
	if (fabs(Ml - Mu) <= 1.0) {
	  zm->q[n]        = (int) (Ml - Mu);
	  zm->shift[n]    = gLl*Ml - gLu*Mu;
          zm->strength[n] = MolZeemanStr(Ju, Mu, Jl, Ml);

	  if ((Ml - Mu) == 1)
	    g_eff += zm->strength[n] * zm->shift[n];

	  norm[zm->q[n]+1] += zm->strength[n];
          n++;
	}
      }
    }
    for (n = 0;  n < zm->Ncomponent;  n++)
      zm->strength[n] /= norm[zm->q[n]+1];
    mrt->g_Lande_eff = g_eff / norm[2];
  }

  vacuum_to_air(1, &(mrt->lambda0), &lambda_air);
  sprintf(messageStr, " -- %2s line at %9.4f nm has %3d "
	  "Zeeman components, gL_eff = %7.4f\n", mrt->molecule->ID,
	  lambda_air, zm->Ncomponent, mrt->g_Lande_eff);
  Error(MESSAGE, routineName, messageStr);

  return zm;
}
/* ------- end ------------------------- MolZeeman.c ---------------- */
