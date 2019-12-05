/* ------- file: -------------------------- h2collisions.c ----------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Sat Sep 19 15:27:09 2009 --

       --------------------------                      ----------RH-- */

/* --- Calculate collisional excitation rates for the vibrational
       transitions in the CO ground electronic state X^1\Sigma^+.
       Calculated are the collisional coefficients Omega_ul [m^3 s^-1],
       where C_ul = Omega_ul * n_pert, and n_pert is the number density
       of perturbers [m^-3].

  See: Ayres, T.R., Wiedemann, G., 1989, ApJ 338, 1033
       --                                              -------------- */


#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "error.h"
#include "constant.h"

/* --- Collisional constants in Landau - Teller format -- ----------- */

#define A_H   3.0
#define B_H  18.1

#define A_HE 87.0
#define B_HE 19.1

#define A_H2 64.0
#define B_H2 19.1

/* --- Harmonic oscillator constant for H2 --          -------------- */

#define OMEGA0 2.103E+05


/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern char messageStr[];


/* ------- begin -------------------------- H2collisions.c ---------- */

void H2collisions(struct Molecule *molecule)
{
  const char routineName[] = "H2collisions";
  register int k;

  double hcomega_k, *beta, C_0;

  if (!strstr(molecule->ID, "H2")) {
    sprintf(messageStr, "Molecule is not H2: %s\n", molecule->ID);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }

  C_0 = KBOLTZMANN / ATM_TO_PA;

  hcomega_k = HPLANCK * CLIGHT * OMEGA0 / KBOLTZMANN;
  beta = (double *) malloc(atmos.Nspace * sizeof(double));
  for (k = 0;  k < atmos.Nspace;  k++)
    beta[k] = hcomega_k / atmos.T[k];

  /* --- For now only calculate a v-independent collisional
         de-excitation rate C_ul. The upward rate is then given by:

           C_lu = nv^*_u / nv^*_l * C_ul

	 --                                            -------------- */

  molecule->C_ul = (double *) malloc(atmos.Nspace * sizeof(double));

  for (k = 0;  k < atmos.Nspace;  k++) {
    if (molecule->n[k]) {

      /* --- Neutral hydrogen contribution --          -------------- */

      molecule->C_ul[k] = C_0 * atmos.T[k] *
	exp(B_H - A_H * pow(atmos.T[k], -0.33333333)) /
	(1.0 - exp(-beta[k])) * atmos.H->n[0][k];

      /* --- Neutral helium contribution
	     (assuming all helium is neutral) --       -------------- */

      molecule->C_ul[k] += C_0 * atmos.T[k] *
	exp(B_HE - A_HE * pow(atmos.T[k], -0.33333333)) /
	(1.0 - exp(-beta[k])) * atmos.elements[1].abund * atmos.nHtot[k];

      /* --- H2 contribution --                        -------------- */

      molecule->C_ul[k] += C_0 * atmos.T[k] *
	exp(B_H2 - A_H2 * pow(atmos.T[k], -0.33333333)) /
	(1.0 - exp(-beta[k])) * atmos.H2->n[k];

      /* --- Electronic contribution --                -------------- */

      molecule->C_ul[k] += (1.4E-09 * CUBE(CM_TO_M)) / sqrt(beta[k]) *
	((1.0 + beta[k]) + 19.0*exp(-3.22*beta[k]) *
	 (1.0 + 4.22*beta[k])) * atmos.ne[k];
    }
  }
  free(beta);
}
/* ------- end ---------------------------- H2collisions.c ---------- */
