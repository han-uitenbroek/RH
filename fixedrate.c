/* ------- file: -------------------------- fixedrate.c -------------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Tue Jan 25 11:51:32 2000 --

       --------------------------                      ----------RH-- */

/* --- Evaluate fixed radiative rates assuming that the radiation field
       is given by the Planck function at the radiation temperature,
       that the radiation filed is constant over the line profile of a
       bound-bound transition, and that the ionization cross section of
       bound-free transitions falls off with a \nu^-3 dependence.

       There are three options for the run of the radiation temperature,
       which is set separately for each fixed transition ft:

       TRAD_ATMOSPHERIC    -- Radiation temperature Trad will follow
                              the kinetic temperature Te throughout the
                              atmosphere.

       TRAD_PHOTOSPHERIC   -- Radiation temperature Trad is equal to kinetic
                              temperature Te in deep photosphere till Te
                              reaches the specified ft->Trad for that
                              transition. It is constant at ft->Trad
                              further outward.

       TRAD_CHROMOSPHERIC  -- Radiation temperature follows the kinetic
                              temperature from deep photosphere through
                              temperature minimum, and then is constant at
                              ft->Trad when Te > ft->Trad in the chromosphere.


      Convention: \Gamma_ij = Gamma[i][j] represents the
                   transition j --> i
       --                                              -------------- */ 


#include <math.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "constant.h"
#include "statistics.h"

#define N_MAX_INTEGRAL  999
#define DELTA_MIN       1.0E-05


/* --- Function prototypes --                          -------------- */

double E1(double x);
double sumE1(double xrad);
double sumE1_xe(double xe, double xrad);


/* --- Global variables --                             -------------- */

extern enum Topology topology;
extern Atmosphere atmos;


/* ------- begin -------------------------- FixedRate.c ------------- */

void FixedRate(Atom *atom)
{
  register int kf, k;

  int    Nrepeat, Ndepth, i, j, ij, ji;
  double Te, Trad, deltaT, gij, e_mc, xe, xrad, C0, hc_kla, exp_hckTla;
  FixedTransition *ft = atom->ft;

  getCPU(3, TIME_START, NULL);

  e_mc = (Q_ELECTRON / EPSILON_0) * (Q_ELECTRON / (CLIGHT*M_ELECTRON));
  if (topology == TWO_D_PLANE) {
    Nrepeat = atmos.N[0];
    Ndepth  = atmos.N[1];
  } else if (topology == THREE_D_PLANE) {
    Nrepeat = atmos.N[0] * atmos.N[1];
    Ndepth  = atmos.N[2];
  } else {
    Nrepeat = 1;
    Ndepth  = atmos.Nspace;
  }
  /* --- Go through the fixed transitions --           -------------- */

  for (kf = 0;  kf < atom->Nfixed;  kf++) {
    ft = atom->ft + kf;

    hc_kla = (HPLANCK * CLIGHT) / (KBOLTZMANN * NM_TO_M * ft->lambda0);
    i = ft->i;
    j = ft->j;
    ij = i*atom->Nlevel + j;
    ji = j*atom->Nlevel + i;

    switch (ft->type) {
    case FIXED_LINE:
      C0  = 2.0*PI * e_mc * ft->strength / SQ(NM_TO_M * ft->lambda0);
      gij = atom->g[i] / atom->g[j];
      break;
    case FIXED_CONTINUUM:
      gij = 0.0;
      C0 = 8.0*PI * ft->strength * CLIGHT / CUBE(NM_TO_M * ft->lambda0);
    }

    for (k = 0;  k < atmos.Nspace;  k++) {
      Te = atmos.T[k];
      if (k >= Nrepeat * (Ndepth - 1))
	deltaT = -1.0;
      else
	deltaT = Te - atmos.T[k + Nrepeat];
      
      /* --- Determine radiation temperature --      -------------- */
      
      Trad = Te;
      switch (ft->option) {
      case TRAD_ATMOSPHERIC:
	break;
      case TRAD_PHOTOSPHERIC:
	if (deltaT >= 0.0  ||  (deltaT < 0.0 && ft->Trad > Te))
	  Trad = ft->Trad;
	break;
      case TRAD_CHROMOSPHERIC:
	if (deltaT >= 0.0  &&  ft->Trad < Te)
	  Trad = ft->Trad;
	break;
      }
      /* --- Evaluate the rates and store in collisional matrix -- -- */      
      
      xrad = hc_kla / Trad;
      switch (ft->type) {
      case FIXED_LINE:
	exp_hckTla = exp(-xrad);
	atom->C[ji][k] += C0 * exp_hckTla / (1.0 - exp_hckTla);
	atom->C[ij][k] += gij * C0 / (1.0 - exp_hckTla);
	break;
      case FIXED_CONTINUUM:
	xe = hc_kla / Te;
	atom->C[ji][k] += C0 * sumE1(xrad);
	atom->C[ij][k] += atom->nstar[i][k] / atom->nstar[j][k] * C0 * 
	  sumE1_xe(xe, xrad);
        break;
      }
    }
  }
  getCPU(3, TIME_POLL, "Fixed Rates");
}
/* ------- end ---------------------------- FixedRate.c ------------- */

/* ------- begin -------------------------- sumE1.c ----------------- */

double sumE1(double xrad)
{
  register int n;

  double sum, dsum;

  /*                     oo
   *  Evaluate the sum  Sum E1(n*xrad)
   *                    n=1
   */

  if ((sum = E1(xrad)) == 0.0) return sum;

  for (n = 2;  n < N_MAX_INTEGRAL;  n++) {
    dsum = E1(n * xrad);
    sum += dsum;
    if (dsum / sum <= DELTA_MIN) break;
  }
  return sum;
}
/* ------- end ---------------------------- sumE1.c ----------------- */

/* ------- begin -------------------------- sumE1_xe.c -------------- */

double sumE1_xe(double xe, double xrad)
{
  register int n;

  double sum, dsum;

  /*                     oo
   *  Evaluate the sum  Sum E1(xe + n*xrad)
   *                    n=0
   */

  if ((sum = E1(xe) + E1(xe + xrad)) == 0.0) return sum;

  for (n = 2;  n < N_MAX_INTEGRAL;  n++) {
    dsum = E1(xe + n * xrad);
    sum += dsum;
    if (dsum / sum <= DELTA_MIN) break;
  }
  return sum;
}
/* ------- end ---------------------------- sumE1_xe.c -------------- */
