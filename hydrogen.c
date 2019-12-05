/* ------- file: -------------------------- hydrogen.c --------------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Thu Jan 26 15:47:01 2012 --

       --------------------------                      ----------RH-- */

/* --- Computes hydrogen related bound-free and free-free opacities,
       and continuum opacities.

       Currently includes:

    --  Neutral Hydrogen bound-free and free-free.
    --  H^- bound-free and free-free.
    --  H2^- free-free.
    --  H2^+ free-free.
    --  Rayleigh scattering by molecular H2.

  Note: For Rayleigh scattering by neutral hydrogen the general
        routine rayleigh.c is used (called from Background).


       Global variables:
        atmos -- Atmos structure for atmospheric data.
        atom  -- Atom structure for current active atomic model.

       Input:
        lambda -- Wavelength [nm] for which opacity and emissivity 
                  are to be calculated.
             H -- Pointer to Atom structure H containing atomic data
                  for hydrogen in atmosphere.

       Output:
        chi[Nspace] -- Array for opacities [m^2].
        eta[Nspace] -- Array for emissivities [J s^-1 Hz^-1 sr^-1].

       --                                              -------------- */
 
#include <ctype.h>
#include <stdlib.h>
#include <math.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "constant.h"
#include "background.h"
#include "error.h"


/* --- Function prototypes --                          -------------- */

double bilinear(int Ncol, int Nrow, double *f, double x, double y);


/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern char messageStr[];


/* ------- begin -------------------------- distribute_nH.c --------- */

void distribute_nH()
{
  const char routineName[] = "distribute_nH";
  register int k, i;

  char    config[4], *ptr;
  int    *quantumNo, Nspace = atmos.Nspace, iq;
  double *g_total, g_i;

  /* --- Redistribute the hydrogen levels for use in the background
         opacity package.

         Use the following conventions:

         1 -- To save memory:
              - let atmos.np point to atmos.H.n[atmos.H.Nlevel-1]
              - let nH2 point to atmos.molecules[0].n;

         2 -- atmos.nHtot represents the total number of hydrogen
              atoms in neutral atoms, protons, H-, and H2 and H2+
              molecules. So in general atmos.nHtot != atmos.H.ntotal.

         3 -- In case of LTE (set by atmos.H_LTE, see keyword.input)
              let atmos.H.n point to the LTE equivalents atmos.H.nstar.

     --                                                -------------- */

  if (atmos.H_LTE) {
    atmos.H->NLTEpops = FALSE;
    Error(MESSAGE, routineName,
	  "\nUsing LTE hydrogen populations for background opacities\n\n");

    /* --- To save memory space let atmos.H->n point to LTE populations
           atmos.H->nstar --                           -------------- */

    if (!atmos.H->active) atmos.H->n = atmos.H->nstar;

  } else {
    atmos.H->NLTEpops = TRUE;
    if (!atmos.H->active)
      atmos.H->n = matrix_double(atmos.H->Nlevel, atmos.Nspace);

    /* --- Find principal quantum number from label --   ------------ */
    
    quantumNo = (int *) malloc((atmos.H->Nlevel - 1) * sizeof(int));
    g_total   = (double *) calloc(atmos.NHydr - 1, sizeof(double));
    
    for (i = 0;  i < atmos.H->Nlevel-1;  i++) {
      sscanf(atmos.H->label[i], "H I %s", config);
      ptr = config;  while (isdigit(*ptr)) ptr++;  *ptr = ' ';
      sscanf(config, "%d", &quantumNo[i]);

      if (quantumNo[i] < atmos.NHydr)
	g_total[quantumNo[i] - 1] += atmos.H->g[i];
    }
    /* --- Now redistribute atmosphere's Hydrogen populations over
           the levels of atom H --                     -------------- */
    
    for (i = 0;  i < atmos.H->Nlevel-1;  i++) {
      if (quantumNo[i] < atmos.NHydr) {
	iq  = quantumNo[i] - 1;
	g_i = atmos.H->g[i] / g_total[iq];
	for (k = 0;  k < Nspace;  k++)
	  atmos.H->n[i][k] = g_i * atmos.nH[iq][k];
      } else {
        sprintf(messageStr, "Too many hydrogen levels (level n = %d)\n"
		" Background opacity additional levels set to zero%s",
		quantumNo[i], (i == atmos.H->Nlevel-2) ? "\n\n" : "");
	Error(WARNING, routineName, messageStr);
      }
    }    
    free(quantumNo);  free(g_total);

    /* --- The protons come last --                    -------------- */

    for (k = 0;  k < Nspace;  k++)
      atmos.H->n[atmos.H->Nlevel-1][k] = atmos.nH[atmos.NHydr-1][k];
  }
  /* --- Free memory for atmospheric populations --    -------------- */

  freeMatrix((void **) atmos.nH);
}
/* ------- end ---------------------------- distribute_nH.c --------- */

/* ------- begin -------------------------- Hydrogen_bf.c ----------- */

bool_t Hydrogen_bf(double lambda, double *chi, double *eta)
{
  /* --- Hydrogen bound-free opacity

    See: Mihalas (1978) p. 99 --                       -------------- */

  register int  k, kr;

  bool_t  opaque;
  int     i;
  double  lambdaEdge, sigma, sigma0, g_bf, twohnu3_c2, twohc, gijk,
    hc_k, hc_kla, *npstar, expla, n_eff, *np;
  AtomicContinuum *continuum;

  opaque = FALSE;
  for (k = 0;  k < atmos.Nspace;  k++) {
    chi[k] = 0.0;
    eta[k] = 0.0;
  }
  if (atmos.H->active) return opaque;

  twohc  = (2.0 * HPLANCK * CLIGHT) / CUBE(NM_TO_M);
  hc_k   = (HPLANCK * CLIGHT) / (KBOLTZMANN * NM_TO_M);
  sigma0 = 32.0/(3.0*sqrt(3.0)) * SQ(Q_ELECTRON)/(4.0*PI*EPSILON_0) / 
    (M_ELECTRON * CLIGHT) * HPLANCK/(2.0*E_RYDBERG);

  npstar = atmos.H->nstar[atmos.H->Nlevel - 1];

  for (kr = 0;  kr < atmos.H->Ncont;  kr++) {
    continuum = atmos.H->continuum + kr;
    lambdaEdge = continuum->lambda0;
    i = continuum->i;
    
    if (lambda <= lambdaEdge  &&  lambda >= continuum->lambda[0]) {
      opaque = TRUE;

      /* --- Find the principal quantum number of level i -- -------- */

      n_eff = sqrt(E_RYDBERG /
		   (atmos.H->E[continuum->j] - atmos.H->E[continuum->i]));

      g_bf  = Gaunt_bf(lambda, n_eff, atmos.H->stage[i] + 1);
      sigma = sigma0 * n_eff * g_bf * CUBE(lambda/lambdaEdge);
      hc_kla     = hc_k / lambda;
      twohnu3_c2 = twohc / CUBE(lambda);
      
      np = atmos.H->n[atmos.H->Nlevel-1];
      for (k = 0;  k < atmos.Nspace;  k++) {
	expla   = exp(-hc_kla/atmos.T[k]);
	gijk    = atmos.H->nstar[i][k]/npstar[k] * expla;
	chi[k] += sigma * (1.0 - expla) * atmos.H->n[i][k];
	eta[k] += twohnu3_c2 * gijk * sigma * np[k];
      }
    }
  }
  return opaque;
}
/* ------- end ---------------------------- Hydrogen_bf.c ----------- */

/* ------- begin -------------------------- Hydrogen_ff.c ----------- */

void Hydrogen_ff(double lambda, double *chi)
{
  /* --- Hydrogen free-free opacity

    See: Mihalas (1978) p. 101
        --                                             -------------- */

  register int k;

  int    Nspace = atmos.Nspace;
  double hc_kla, C0, sigma, g_ff, stim, nu3, *np;

  C0     = SQ(Q_ELECTRON)/(4.0*PI*EPSILON_0) / sqrt(M_ELECTRON);
  sigma  = 4.0/3.0 * sqrt(2.0*PI/(3.0 * KBOLTZMANN)) * CUBE(C0) /
    (HPLANCK * CLIGHT);
  nu3    = CUBE((lambda * NM_TO_M) / CLIGHT);
  hc_kla = (HPLANCK * CLIGHT) / (KBOLTZMANN * NM_TO_M * lambda);

  np = atmos.H->n[atmos.H->Nlevel-1];
  for (k = 0;  k < Nspace;  k++) {
    stim   = 1.0 - exp(-hc_kla/atmos.T[k]);
    g_ff   = Gaunt_ff(lambda, 1, atmos.T[k]);
    chi[k] = sigma / sqrt(atmos.T[k]) * nu3 * atmos.ne[k] *
      np[k] * stim * g_ff;
  }
}
/* ------- end ---------------------------- Hydrogen_ff.c ----------- */

/* ------- begin -------------------------- Gaunt_bf.c -------------- */

double Gaunt_bf(double lambda, double n_eff, int charge)
{
  /* --- M. J. Seaton (1960), Rep. Prog. Phys. 23, 313 -- ----------- */

  double x, x3, nsqx;

  x    = ((HPLANCK*CLIGHT)/(lambda * NM_TO_M)) / (E_RYDBERG * SQ(charge));
  x3   = pow(x, 0.33333333);
  nsqx = 1.0 / (SQ(n_eff) * x);

  return 1.0 + 0.1728*x3 * (1.0 - 2.0*nsqx) -
               0.0496*SQ(x3) * (1.0 - (1.0 - nsqx)*0.66666667*nsqx);
}
/* ------- end ---------------------------- Gaunt_bf.c -------------- */

/* ------- begin -------------------------- Gaunt_ff.c -------------- */

double Gaunt_ff(double lambda, int charge, double T)
{
  /* --- M. J. Seaton (1960), Rep. Prog. Phys. 23, 313

   Note: There is a problem with this expansion at higher temperatures
         (T > 3.0E4 and longer wavelengths (lambda > 2000 nm). Set to
         1.0 when the value goes below 1.0 --          -------------- */

  double x, x3, y, gIII;

  x  = ((HPLANCK * CLIGHT)/(lambda * NM_TO_M)) / (E_RYDBERG * SQ(charge));
  x3 = pow(x, 0.33333333);
  y  = (2.0 * lambda * NM_TO_M * KBOLTZMANN*T) / (HPLANCK*CLIGHT);

  gIII = 1.0 + 0.1728*x3 * (1.0 + y) -
               0.0496*SQ(x3) * (1.0 + (1.0 + y)*0.33333333*y);
  return (gIII > 1.0) ? gIII : 1.0;
}
/* ------- end ---------------------------- Gaunt_ff.c -------------- */

/* ------- begin -------------------------- Hminus_bf.c ------------- */

#define NBF            34

bool_t Hminus_bf(double lambda, double *chi, double *eta)
{
  register int k;

  /* --- H-minus Bound-Free coefficients (in units of 1.0E-21 m^2).

   From: S. Geltman (1962), ApJ 136, 935-945
   Also: Mihalas (1978), p. 102 --                     -------------- */

  static double lambdaBF[NBF] = {
    0.0, 50.0, 100.0, 150.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0,
    500.0, 550.0, 600.0, 650.0, 700.0, 750.0, 800.0, 850.0, 900.0, 950.0,
    1000.0, 1050.0, 1100.0, 1150.0, 1200.0, 1250.0, 1300.0, 1350.0, 
    1400.0, 1450.0, 1500.0, 1550.0, 1600.0, 1641.9};

  static double alphaBF[NBF] = {
    0.0,  0.15, 0.33, 0.57, 0.85, 1.17, 1.52, 1.89, 2.23, 2.55, 2.84,
    3.11, 3.35, 3.56, 3.71, 3.83, 3.92, 3.95, 3.93, 3.85, 3.73, 3.58,
    3.38, 3.14, 2.85, 2.54, 2.20, 1.83, 1.46, 1.06, 0.71, 0.40, 0.17, 0.0};

  bool_t hunt;
  int    Nspace = atmos.Nspace;
  double hc_kla, stimEmis, twohnu3_c2, alpha_bf;

  if ((lambda <= lambdaBF[0]) || (lambda >= lambdaBF[NBF-1]))
    return FALSE;

  splineCoef(NBF, lambdaBF, alphaBF);
  splineEval(1, &lambda, &alpha_bf, hunt=FALSE);
  alpha_bf *= 1.0E-21;

  hc_kla     = (HPLANCK * CLIGHT) / (KBOLTZMANN * NM_TO_M * lambda);
  twohnu3_c2 = (2.0 * HPLANCK * CLIGHT) / CUBE(NM_TO_M * lambda);

  for (k = 0;  k < Nspace;  k++) {
    stimEmis = exp(-hc_kla/atmos.T[k]);
    chi[k]   = atmos.nHmin[k] * (1.0 - stimEmis) * alpha_bf;
    eta[k]   = atmos.nHmin[k] * twohnu3_c2 * stimEmis * alpha_bf;
  }

  return TRUE;
}
/* ------- end ---------------------------- Hminus_bf.c ------------- */

/* ------- begin -------------------------- Hminus_ff.c ------------- */

#define NFF            17
#define NTHETA         16

bool_t Hminus_ff(double lambda, double *chi)
{
  register int  k;

  static bool_t  initialize=TRUE;
  static int     index;
  static double *theta_index;

  /* --- H-minus Free-Free coefficients (in units of 1.0E-29 m^5/J)

   From: J. L. Stilley and J. Callaway (1970), ApJ 160, 245-260
   Also: D. Mihalas (1978), p. 102
         R. Mathisen (1984), Master's thesis, Inst. Theor.
          Astroph., University of Oslo. p. 17

         When called the first time (or when initialize==TRUE) the
         fractional indices for atmospheric temperatures into the
         theta table are stored in theta_index. This memory can be
         freed by calling the routine with lambda==0.0
         --                                            -------------- */

  static double lambdaFF[NFF] = {
    0.0, 303.8, 455.6, 506.3, 569.5, 650.9, 759.4, 911.3, 1013.0, 1139.0,
    1302.0, 1519.0, 1823.0, 2278.0, 3038.0, 4556.0, 9113.0};

  /* --- theta = 5040.0/T --                           -------------- */

  static double thetaFF[NTHETA] = {
    0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2,
    1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0};

  static double kappaFF[NFF * NTHETA] = {
    /* --- lambda =    0.0 [nm] --                     -------------- */
     0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
    /* --- lambda =  303.8 [nm] --                     -------------- */
     3.44e-02, 4.18e-02, 4.91e-02, 5.65e-02, 6.39e-02, 7.13e-02,
     7.87e-02, 8.62e-02, 9.36e-02, 1.01e-01, 1.08e-01, 1.16e-01,
     1.23e-01, 1.30e-01, 1.38e-01, 1.45e-01,
    /* --- lambda =  455.6 [nm] --                     -------------- */
     7.80e-02, 9.41e-02, 1.10e-01, 1.25e-01, 1.40e-01, 1.56e-01,
     1.71e-01, 1.86e-01, 2.01e-01, 2.16e-01, 2.31e-01, 2.45e-01,
     2.60e-01, 2.75e-01, 2.89e-01, 3.03e-01,
    /* --- lambda =  506.3 [nm] --                     -------------- */
     9.59e-02, 1.16e-01, 1.35e-01, 1.53e-01, 1.72e-01, 1.90e-01,
     2.08e-01, 2.25e-01, 2.43e-01, 2.61e-01, 2.78e-01, 2.96e-01,
     3.13e-01, 3.30e-01, 3.47e-01, 3.64e-01,
    /* --- lambda =  569.5 [nm] --                     -------------- */
     1.21e-01, 1.45e-01, 1.69e-01, 1.92e-01, 2.14e-01, 2.36e-01,
     2.58e-01, 2.80e-01, 3.01e-01, 3.22e-01, 3.43e-01, 3.64e-01,
     3.85e-01, 4.06e-01, 4.26e-01, 4.46e-01,
    /* --- lambda =  650.9 [nm] --                     -------------- */
     1.56e-01, 1.88e-01, 2.18e-01, 2.47e-01, 2.76e-01, 3.03e-01,
     3.31e-01, 3.57e-01, 3.84e-01, 4.10e-01, 4.36e-01, 4.62e-01,
     4.87e-01, 5.12e-01, 5.37e-01, 5.62e-01,
    /* --- lambda =  759.4 [nm] --                     -------------- */
     2.10e-01, 2.53e-01, 2.93e-01, 3.32e-01, 3.69e-01, 4.06e-01,
     4.41e-01, 4.75e-01, 5.09e-01, 5.43e-01, 5.76e-01, 6.08e-01,
     6.40e-01, 6.72e-01, 7.03e-01, 7.34e-01,
    /* --- lambda =  911.3 [nm] --                     -------------- */
     2.98e-01, 3.59e-01, 4.16e-01, 4.70e-01, 5.22e-01, 5.73e-01,
     6.21e-01, 6.68e-01, 7.15e-01, 7.60e-01, 8.04e-01, 8.47e-01,
     8.90e-01, 9.32e-01, 9.73e-01, 1.01e+00,
    /* --- lambda = 1013.0 [nm] --                     -------------- */
     3.65e-01, 4.39e-01, 5.09e-01, 5.75e-01, 6.39e-01, 7.00e-01,
     7.58e-01, 8.15e-01, 8.71e-01, 9.25e-01, 9.77e-01, 1.03e+00,
     1.08e+00, 1.13e+00, 1.18e+00, 1.23e+00,
    /* --- lambda = 1139.0 [nm] --                     -------------- */
     4.58e-01, 5.50e-01, 6.37e-01, 7.21e-01, 8.00e-01, 8.76e-01,
     9.49e-01, 1.02e+00, 1.09e+00, 1.15e+00, 1.22e+00, 1.28e+00,
     1.34e+00, 1.40e+00, 1.46e+00, 1.52e+00,
    /* --- lambda = 1302.0 [nm] --                     -------------- */
     5.92e-01, 7.11e-01, 8.24e-01, 9.31e-01, 1.03e+00, 1.13e+00,
     1.23e+00, 1.32e+00, 1.40e+00, 1.49e+00, 1.57e+00, 1.65e+00,
     1.73e+00, 1.80e+00, 1.88e+00, 1.95e+00,
    /* --- lambda = 1519.0 [nm] --                     -------------- */
     7.98e-01, 9.58e-01, 1.11e+00, 1.25e+00, 1.39e+00, 1.52e+00,
     1.65e+00, 1.77e+00, 1.89e+00, 2.00e+00, 2.11e+00, 2.21e+00,
     2.32e+00, 2.42e+00, 2.51e+00, 2.61e+00,
    /* --- lambda = 1823.0 [nm] --                     -------------- */
     1.14e+00, 1.36e+00, 1.58e+00, 1.78e+00, 1.98e+00, 2.17e+00,
     2.34e+00, 2.52e+00, 2.68e+00, 2.84e+00, 3.00e+00, 3.15e+00,
     3.29e+00, 3.43e+00, 3.57e+00, 3.70e+00,
    /* --- lambda = 2278.0 [nm] --                     -------------- */
     1.77e+00, 2.11e+00, 2.44e+00, 2.75e+00, 3.05e+00, 3.34e+00,
     3.62e+00, 3.89e+00, 4.14e+00, 4.39e+00, 4.63e+00, 4.86e+00,
     5.08e+00, 5.30e+00, 5.51e+00, 5.71e+00,
    /* --- lambda = 3038.0 [nm] --                     -------------- */
     3.10e+00, 3.71e+00, 4.29e+00, 4.84e+00, 5.37e+00, 5.87e+00,
     6.36e+00, 6.83e+00, 7.28e+00, 7.72e+00, 8.14e+00, 8.55e+00,
     8.95e+00, 9.33e+00, 9.71e+00, 1.01e+01,
    /* --- lambda = 4556.0 [nm] --                     -------------- */
     6.92e+00, 8.27e+00, 9.56e+00, 1.08e+01, 1.19e+01, 1.31e+01,
     1.42e+01, 1.52e+01, 1.62e+01, 1.72e+01, 1.82e+01, 1.91e+01,
     2.00e+01, 2.09e+01, 2.17e+01, 2.25e+01,
    /* --- lambda = 9113.0 [nm] --                     -------------- */
     2.75e+01, 3.29e+01, 3.80e+01, 4.28e+01, 4.75e+01, 5.19e+01,
     5.62e+01, 6.04e+01, 6.45e+01, 6.84e+01, 7.23e+01, 7.60e+01,
     7.97e+01, 8.32e+01, 8.67e+01, 9.01e+01
  };
  
  int     Nspace = atmos.Nspace;
  double  theta, pe, lambda_index, kappa;

  if (lambda == 0.0) {

    /* --- When called with zero wavelength free memory for fractional
           indices --                                  -------------- */

    if (theta_index) free(theta_index);
    initialize = TRUE;
    return FALSE;
  }
  /* --- Use long-wavelength expansion if wavelength beyond 9113 nm - */

  if (lambda >= lambdaFF[NFF-1])
    return Hminus_ff_long(lambda, chi);

  if (initialize) {

    /* --- Store the fractional indices of temperature only the
           first time around --                        -------------- */

    theta_index = (double *) malloc(Nspace * sizeof(double));
    for (k = 0;  k < Nspace;  k++) {
      theta = THETA0 / atmos.T[k];
      if (theta <= thetaFF[0])
        theta_index[k] = 0;
      else if (theta >= thetaFF[NTHETA-1])
	theta_index[k] = NTHETA - 1;
      else {
	Hunt(NTHETA, thetaFF, theta, &index);
	theta_index[k] = (double) index +
	  (theta - thetaFF[index]) / (thetaFF[index+1] - thetaFF[index]);
      }
    }
    initialize = FALSE;
  }

  Hunt(NFF, lambdaFF, lambda, &index);
  lambda_index = (double) index +
    (lambda - lambdaFF[index]) / (lambdaFF[index+1] - lambdaFF[index]);

  for (k = 0;  k < Nspace;  k++) {
    pe     = atmos.ne[k] * KBOLTZMANN * atmos.T[k];
    kappa  = bilinear(NTHETA, NFF, kappaFF,
		      theta_index[k], lambda_index);
    chi[k] = (atmos.H->n[0][k] * 1.0E-29) * pe * kappa;
  }
  return TRUE;
}
/* ------- end ---------------------------- Hminus_ff.c ------------- */

#define NJOHN  6

/* ------- begin -------------------------- Hminus_ff_long.c -------- */

bool_t Hminus_ff_long(double lambda, double *chi)
{
  register int k, n;

  /* --- H-minus Free-Free opacity. Parametrization for long wavelengths
         as given by T. L. John (1988), A&A 193, 189-192 (see table 3a).
	 His results are based on calculations by K. L. Bell and
	 K. A. Berrington (1987), J. Phys. B 20, 801-806. -- -------- */

  static double A[NJOHN] = {     0.000,  2483.346, -3449.889,  2200.040,
		              -696.271,    88.283 };
  static double B[NJOHN] = {     0.000,   285.827, -1158.382,  2427.719,
			     -1841.400,   444.517 };
  static double C[NJOHN] = {     0.000, -2054.291,  8746.523,-13651.105,
                              8624.970, -1863.864 };
  static double D[NJOHN] = {     0.000,  2827.776,-11485.632, 16755.524,
                            -10051.530,  2095.288 };
  static double E[NJOHN] = {     0.000, -1341.537,  5303.609, -7510.494,
                              4400.067,  -901.788 };
  static double F[NJOHN] = {     0.000,   208.952,  -812.939,  1132.738,
                              -655.020,   132.985 };

  double  Clambda[NJOHN], lambda_mu, lambda_inv, sqrt_theta, theta_n,
          Ck = (KBOLTZMANN * THETA0 * 1.0E-32);
  
  /* --- First evaluate the wavelength dependent coefficients -- ---- */

  lambda_mu  = lambda / MICRON_TO_NM;
  lambda_inv = 1.0 / lambda_mu;
  for (n = 1;  n < NJOHN;  n++) {
    Clambda[n] = SQ(lambda_mu)*A[n] + B[n] + lambda_inv*(C[n] +
                 lambda_inv*(D[n] + lambda_inv*(E[n] + lambda_inv*F[n])));
  }
  /* --- Then spatial dependence --                      ------------ */

  for (k = 0;  k < atmos.Nspace;  k++) {
    chi[k]     = 0.0;
    theta_n    = 1.0;
    sqrt_theta = sqrt(THETA0 / atmos.T[k]);
    for (n = 1;  n < NJOHN;  n++) {
      theta_n *= sqrt_theta;
      chi[k]  += theta_n * Clambda[n];
    }
    chi[k] *= atmos.H->n[0][k] * (atmos.ne[k] * Ck);
  }
  return TRUE;
}
/* ------- end ---------------------------- Hminus_ff_long.c -------- */

/* ------- begin -------------------------- H2minus_ff.c ------------ */

#define	NFF_H2	    19
#define	NTHETA_H2    8

bool_t H2minus_ff(double lambda, double *chi) {

  register int  k;

  static  bool_t initialize=TRUE;
  static   int   index;
  static double *theta_index;

  /* --- H2-minus Free-Free absorption coefficients (in units of
         10E-29 m^5/J). Stimulated emission is included.

   From: Bell, K. L., (1980) J. Phys. B13, 1859.
   Also: R. Mathisen (1984), Master's thesis, Inst. Theor.
          Astroph., University of Oslo, p. 18

         When called the first time (or when initialize==TRUE) the
         fractional indices for atmospheric temperatures into the
         theta table are stored in theta_index. This memory can be
         freed by calling the routine with lambda==0.0
         --                                            -------------- */


  static double lambdaFF[NFF_H2] = {
       0.0,   350.5,   414.2,   506.3,   569.6,   650.9,   759.4,   911.3,
    1139.1,  1518.8,  1822.6,  2278.3,  3037.7,  3645.2,  4556.5,  6075.3,
    9113.0, 11391.3, 15188.3};

  static double thetaFF[NTHETA_H2] = {
    0.5, 0.8, 1.0, 1.2, 1.6, 2.0, 2.8, 3.6};

  static double kappaFF[NFF_H2 * NTHETA_H2] = {
    /* --- lambda =     0.0 [nm] --                    -------------- */
     0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
     0.00e+00, 0.00e+00, 0.00e+00,
    /* --- lambda =   350.5 [nm] --                    -------------- */
     4.17e-02, 6.10e-02, 7.34e-02, 8.59e-02, 1.11e-01,
     1.37e-01, 1.87e-01, 2.40e-01,
    /* --- lambda =   414.2 [nm] --                    -------------- */
     5.84e-02, 8.43e-02, 1.01e-01, 1.17e-01, 1.49e-01,
     1.82e-01, 2.49e-01, 3.16e-01,
    /* --- lambda =   506.3 [nm] --                    -------------- */
     8.70e-02, 1.24e-01, 1.46e-01, 1.67e-01, 2.10e-01,
     2.53e-01, 3.39e-01, 4.27e-01,
    /* --- lambda =   569.6 [nm] --                    -------------- */
     1.10e-01, 1.54e-01, 1.80e-01, 2.06e-01, 2.55e-01,
     3.05e-01, 4.06e-01, 5.07e-01,
    /* --- lambda =   650.9 [nm] --                    -------------- */
     1.43e-01, 1.98e-01, 2.30e-01, 2.59e-01, 3.17e-01,
     3.75e-01, 4.92e-01, 6.09e-01,
    /* --- lambda =   759.4 [nm] --                    -------------- */
     1.92e-01, 2.64e-01, 3.03e-01, 3.39e-01, 4.08e-01,
     4.76e-01, 6.13e-01, 7.51e-01,
    /* --- lambda =   911.3 [nm] --                    -------------- */
     2.73e-01, 3.71e-01, 4.22e-01, 4.67e-01, 5.52e-01,
     6.33e-01, 7.97e-01, 9.63e-01,
    /* --- lambda =  1139.1 [nm] --                    -------------- */
     4.20e-01, 5.64e-01, 6.35e-01, 6.97e-01, 8.06e-01,
     9.09e-01, 1.11e+00, 1.32e+00,
    /* --- lambda =  1518.8 [nm] --                    -------------- */
     7.36e-01, 9.75e-01, 1.09e+00, 1.18e+00, 1.34e+00, 
     1.48e+00, 1.74e+00, 2.01e+00,
    /* --- lambda =  1822.6 [nm] --                    -------------- */
     1.05e+00, 1.39e+00, 1.54e+00, 1.66e+00, 1.87e+00,
     2.04e+00, 2.36e+00, 2.68e+00,
    /* --- lambda =  2278.3 [nm] --                    -------------- */
     1.63e+00, 2.14e+00, 2.36e+00, 2.55e+00, 2.84e+00,
     3.07e+00, 3.49e+00, 3.90e+00,
    /* --- lambda =  3037.7 [nm] --                    -------------- */
     2.89e+00, 3.76e+00, 4.14e+00, 4.44e+00, 4.91e+00,
     5.28e+00, 5.90e+00, 6.44e+00,
    /* --- lambda =  3645.2 [nm] --                    -------------- */
     4.15e+00, 5.38e+00, 5.92e+00, 6.35e+00, 6.99e+00,
     7.50e+00, 8.32e+00, 9.02e+00,
    /* --- lambda =  4556.5 [nm] --                    -------------- */
     6.47e+00, 8.37e+00, 9.20e+00, 9.84e+00, 1.08e+01,
     1.16e+01, 1.28e+01, 1.38e+01,
    /* --- lambda =  6075.3 [nm] --                    -------------- */
     1.15e+01,1.48e+01, 1.63e+01, 1.74e+01, 1.91e+01, 
     2.04e+01, 2.24e+01, 2.40e+01,
    /* --- lambda =  9113.0 [nm] --                    -------------- */
     2.58e+01, 3.33e+01, 3.65e+01, 3.90e+01, 4.27e+01, 
     4.54e+01, 4.98e+01, 5.33e+01,
    /* --- lambda = 11391.3 [nm] --                    -------------- */
     4.03e+01, 5.20e+01, 5.70e+01, 6.08e+01, 6.65e+01, 
     7.08e+01, 7.76e+01, 8.30e+01,
    /* --- lambda = 15188.3 [nm] --                    -------------- */
     7.16e+01, 9.23e+01, 1.01e+02, 1.08e+02, 1.18e+02, 
     1.26e+02, 1.38e+02, 1.47e+02
  };

  int     Nspace = atmos.Nspace;
  double  theta, pe, lambda_index, kappa, *nH2;

  if (lambda == 0.0) {

    /* --- When called with zero wavelength free memory for fractional
           indices --                                  -------------- */

    if (theta_index) free(theta_index);
    initialize = TRUE;
    return FALSE;
  }
  if (lambda >= lambdaFF[NFF_H2-1])
    return FALSE;

  if (initialize) {
    theta_index = (double *) malloc(Nspace * sizeof(double));
    for (k = 0;  k < Nspace;  k++) {
      theta = THETA0 / atmos.T[k];
      if (theta <= thetaFF[0])
        theta_index[k] = 0;
      else if (theta >= thetaFF[NTHETA_H2-1])
	theta_index[k] = NTHETA_H2-1;
      else {
	Hunt(NTHETA_H2, thetaFF, theta, &index);
	theta_index[k] = index + (theta - thetaFF[index]) /
	  (thetaFF[index+1] - thetaFF[index]);
      }
    }
    initialize = FALSE;
  }

  Hunt(NFF_H2, lambdaFF, lambda, &index);
  lambda_index = index + (lambda - lambdaFF[index]) /
    (lambdaFF[index+1] - lambdaFF[index]);

  nH2 = atmos.H2->n;
  for (k = 0;  k < Nspace;  k++) {
    if (nH2[k] > 0.0) {
      pe     = atmos.ne[k] * KBOLTZMANN * atmos.T[k];
      kappa  = bilinear(NTHETA_H2, NFF_H2, kappaFF,
			theta_index[k], lambda_index);
      chi[k] = (nH2[k] * 1.0E-29) * pe * kappa;
    } else
      chi[k] = 0.0;
  } 
  return TRUE;
}
/* ------- end ---------------------------- H2minus_ff.c ------------ */

/* ------- begin -------------------------- H2plus_ff.c ------------- */

#define NFF_H2P    15
#define NTEMP_H2P  10

bool_t H2plus_ff(double lambda, double *chi)
{
  register int  k;

  static  bool_t initialize=TRUE;
  static   int   index;
  static double *temp_index;

  /* --- H2+ Free-Free scattering coefficients in units of 
         1.0E-49 m^-1 / (H atom/m^3) / (proton/M^3). Stimulated emission
	 is included. This represents the following interaction:

	   H + H^+ + \nu ---> H + H^+

   From: D. R. Bates (1952), MNRAS 112, 40-44
   Also: R. Mathisen (1984), Master's thesis, Inst. Theor.
          Astroph., University of Oslo, p. 45

         When called the first time (or when initialize==TRUE) the
         fractional indices for atmospheric temperatures into the
         temperature table are stored in temp_index. This memory can be
         freed by calling the routine with lambda==0.0
         --                                            -------------- */


  static double lambdaFF[NFF_H2P] = { 
        0.0,  384.6,  555.6,  833.3, 1111.1, 1428.6, 1666.7,
     2000.0, 2500.0, 2857.1, 3333.3, 4000.0, 5000.0, 6666.7, 10000.0};

  static double tempFF[NTEMP_H2P] = {
    2.5E+03, 3.0E+03, 3.5E+03, 4.0E+03, 5.0E+03, 
    6.0E+03, 7.0E+03, 8.0E+03, 1.0E+04, 1.2E+04};

  static double kappaFF[NFF_H2P * NTEMP_H2P] = {
    /* --- lambda =      0.0 [nm] --                   -------------- */
    0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
    /* --- lambda =    384.6 [nm] --                   -------------- */
    0.46, 0.46, 0.42, 0.39, 0.36, 0.33, 0.32, 0.30, 0.27, 0.25,
    /* --- lambda =    555.6 [nm] --                   -------------- */
    0.70, 0.62, 0.59, 0.56, 0.51, 0.43, 0.41, 0.39, 0.35, 0.34,
    /* --- lambda =    833.3 [nm] --                   -------------- */
    0.92, 0.86, 0.80, 0.76, 0.70, 0.64, 0.59, 0.55, 0.48, 0.43,
    /* --- lambda =   1111.1 [nm] --                   -------------- */
    1.11, 1.04, 0.96, 0.91, 0.82, 0.74, 0.68, 0.62, 0.53, 0.46,
    /* --- lambda =   1428.6 [nm] --                   -------------- */
    1.26, 1.19, 1.09, 1.02, 0.90, 0.80, 0.72, 0.66, 0.55, 0.48,
    /* --- lambda =   1666.7 [nm] --                   -------------- */
    1.37, 1.25, 1.15, 1.07, 0.93, 0.83, 0.74, 0.67, 0.56, 0.49,
    /* --- lambda =   2000.0 [nm] --                   -------------- */
    1.44, 1.32, 1.21, 1.12, 0.97, 0.84, 0.75, 0.67, 0.56, 0.48,
    /* --- lambda =   2500.0 [nm] --                   -------------- */
    1.54, 1.39, 1.26, 1.15, 0.98, 0.85, 0.75, 0.67, 0.55, 0.46,
    /* --- lambda =   2857.1 [nm] --                   -------------- */
    1.58, 1.42, 1.27, 1.16, 0.98, 0.84, 0.74, 0.66, 0.54, 0.45,
    /* --- lambda =   3333.3 [nm] --                   -------------- */
    1.62, 1.43, 1.28, 1.15, 0.97, 0.83, 0.72, 0.64, 0.52, 0.44,
    /* --- lambda =   4000.0 [nm] --                   -------------- */
    1.63, 1.43, 1.27, 1.14, 0.95, 0.80, 0.70, 0.62, 0.50, 0.42,
    /* --- lambda =   5000.0 [nm] --                   -------------- */
    1.62, 1.40, 1.23, 1.10, 0.90, 0.77, 0.66, 0.59, 0.48, 0.39,
    /* --- lambda =   6666.7 [nm] --                   -------------- */
    1.55, 1.33, 1.16, 1.03, 0.84, 0.71, 0.60, 0.53, 0.43, 0.36,
    /* --- lambda =  10000.0 [nm] --                   -------------- */
    1.39, 1.18, 1.02, 0.90, 0.73, 0.60, 0.52, 0.46, 0.37, 0.31
  };

  int     Nspace = atmos.Nspace;
  double  T, lambda_index, kappa, *np;

  if (lambda == 0.0) {

    /* --- When called with zero wavelength free memory for fractional
           indices --                                  -------------- */

    if (temp_index) free(temp_index);
    initialize = TRUE;
    return FALSE;
  }
  if (lambda >= lambdaFF[NFF_H2P-1])
    return FALSE;

  if (initialize) {
    temp_index = (double *) malloc(Nspace * sizeof(double));
    for (k = 0;  k < Nspace;  k++) {
      T = atmos.T[k];
      if (T <= tempFF[0])
        temp_index[k] = 0;
      else if (T >= tempFF[NTEMP_H2P-1])
	temp_index[k] = NTEMP_H2P-1;
      else {
	Hunt(NTEMP_H2P, tempFF, T, &index);
	temp_index[k] = index + (T - tempFF[index]) /
	  (tempFF[index+1] - tempFF[index]);
      }
    }
    initialize = FALSE;
  }

  Hunt(NFF_H2P, lambdaFF, lambda, &index);
  lambda_index = index + (lambda - lambdaFF[index]) /
    (lambdaFF[index+1] - lambdaFF[index]);

  np = atmos.H->n[atmos.H->Nlevel-1];    
  for (k = 0;  k < Nspace;  k++) {
    kappa  = bilinear(NTEMP_H2P, NFF_H2P, kappaFF,
		      temp_index[k], lambda_index);
    chi[k] = (atmos.H->n[0][k] * 1.0E-29) * (np[k] * 1.0E-20) * kappa;
  }

  return TRUE;
}
/* ------- end ---------------------------- H2plus_ff.c ------------- */

/* ------- begin -------------------------- Rayleigh_H2.c ----------- */

#define RAYLEIGH_H2_LIMIT  121.57
#define N_RAYLEIGH_H2      21

bool_t Rayleigh_H2(double lambda, double *scatt)
{
  register int k;

  static double a[3] = {8.779E+01, 1.323E+06, 2.245E+10};

  static double lambdaRH2[N_RAYLEIGH_H2] = {
    121.57, 130.00, 140.00, 150.00, 160.00, 170.00, 185.46,
    186.27, 193.58, 199.05, 230.29, 237.91, 253.56, 275.36,
    296.81, 334.24, 404.77, 407.90, 435.96, 546.23, 632.80 };

  static double sigma[N_RAYLEIGH_H2] = {
    2.35E-06, 1.22E-06, 6.80E-07, 4.24E-07, 2.84E-07, 2.00E-07, 1.25E-07,
    1.22E-07, 1.00E-07, 8.70E-08, 4.29E-08, 3.68E-08, 2.75E-08, 1.89E-08,
    1.36E-08, 8.11E-09, 3.60E-09, 3.48E-09, 2.64E-09, 1.04E-09, 5.69E-10 };

  /* --- Rayleigh scattering by H2 molecules. Cross-section is given
         in in units of Mb, 1.0E-22 m^2.

    See: G. A. Victor and A. Dalgarno (1969), J. Chem. Phys. 50, 2535
         (for lambda <= 632.80 nm), and
         S. P. Tarafdar and M. S. Vardya (1973), MNRAS 163, 261
   Also: R. Mathisen (1984), Master's thesis, Inst. Theor.
         Astroph., University of Oslo, p. 49
         --                                            -------------- */

  bool_t hunt;
  double lambda2, sigma_RH2, *nH2;

  nH2 = atmos.H2->n;

  if (lambda >= RAYLEIGH_H2_LIMIT) {
    if (lambda <= lambdaRH2[N_RAYLEIGH_H2 - 1]) {
      Linear(N_RAYLEIGH_H2, lambdaRH2, sigma,
		   1, &lambda, &sigma_RH2, hunt=FALSE);
    } else {
      lambda2 = 1.0 / SQ(lambda);
      sigma_RH2 = (a[0] + (a[1] + a[2]*lambda2) * lambda2) * SQ(lambda2);
    }
    sigma_RH2 *= MEGABARN_TO_M2;

    for (k = 0;  k < atmos.Nspace;  k++)
      scatt[k] = sigma_RH2 * nH2[k];

    return TRUE;
  } else 
    return FALSE;
}
/* ------- end ---------------------------- Rayleigh_H2.c ----------- */

/* ------- begin -------------------------- bilinear.c -------------- */

double bilinear(int Ncol, int Nrow, double *f, double x, double y)
{
  int    i, j, i1, j1;
  double fx, fy;

  /* --- Bilinear interpolation of the function f on the fractional
         indices x and y --                            -------------- */

  i = (int) x;  fx = x - i;
  if (i == Ncol-1)
    i1 = i;
  else
    i1 = i + 1;
  j = (int) y;  fy = y - j;
  if (j == Nrow-1)
    j1 = j;
  else
    j1 = j + 1;

  return (1.0 - fx)*(1.0 - fy) * f[j*Ncol+i] +
                 fx*(1.0 - fy) * f[j*Ncol+i1] +
                 (1.0 - fx)*fy * f[j1*Ncol+i] +
                         fx*fy * f[j1*Ncol+i1];
}
/* ------- end ---------------------------- bilinear.c -------------- */
