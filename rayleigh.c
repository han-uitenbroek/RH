/* ------- file: -------------------------- rayleigh.c --------------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Wed Jan 26 11:38:53 2000 --

       --------------------------                      ----------RH-- */

#include <stdlib.h>
#include <math.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "constant.h"
#include "background.h"
#include "error.h"

#define  LONG_WAVELENGTH  1.0E6


/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern char messageStr[];


/* ------- begin -------------------------- Rayleigh.c -------------- */

bool_t Rayleigh(double lambda, Atom *atom, double *scatt)
{
  /* --- Rayleigh scattering by transitions from the ground state of
         neutral atoms. Sums scattering crosssections of all bound-bound
         transitions from the groundstate of the atom with lamda_red
         (the longest wavelength covered by the transition) less than
         wavelength lambda.

    See: Mihalas (1978) p. 106 --                      -------------- */

  const char routineName[] = "Rayleigh";
  register int k, kr;

  double lambda_limit, lambda_red, C, sigma_e, fomega, f, lambda2,
         sigma_Rayleigh;
  AtomicLine *line;

  if (atom->stage[0] != 0) {
    sprintf(messageStr, "Lowest level of atom is not from neutral stage: %s",
	    atom->label[0]);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }
  /* --- Find the red limit of the shortest wavelength b-b transition
         from the atom's ground state --               -------------- */

  if (atom->Nline > 0) {
    lambda_limit = LONG_WAVELENGTH;
    for (kr = 0;  kr < atom->Nline;  kr++) {
      line = atom->line + kr;
      if (line->i == 0) {
	lambda_red = line->lambda0 * (1.0 + line->qwing * 
				    atmos.vmicro_char / CLIGHT);
	lambda_limit = MIN(lambda_limit, lambda_red);
      }
    }
  } else
    return FALSE;

  if (lambda > lambda_limit) {
    C = 2*PI * (Q_ELECTRON/EPSILON_0) * (Q_ELECTRON/M_ELECTRON) / CLIGHT;
    sigma_e = 8.0*PI/3.0 * pow(Q_ELECTRON/(sqrt(4.0*PI*EPSILON_0) *
					   (sqrt(M_ELECTRON)*CLIGHT)), 4);

    fomega = 0.0;
    for (kr = 0;  kr < atom->Nline;  kr++) {
      line = atom->line + kr;
      lambda_red = line->lambda0 * (1.0 + line->qwing * 
				  atmos.vmicro_char / CLIGHT);
      if (line->i == 0  &&  lambda > lambda_red) {
	lambda2 = 1.0 / (SQ(lambda / line->lambda0) - 1.0);
        f = line->Aji * (atom->g[line->j] / atom->g[0]) *
	  SQ(line->lambda0*NM_TO_M) / C;
	fomega += f * SQ(lambda2);
      }
    }
    sigma_Rayleigh = sigma_e * fomega;

    for (k = 0;  k < atmos.Nspace;  k++)
      scatt[k] = sigma_Rayleigh * atom->n[0][k];

    return TRUE;
  } else
    return FALSE;
}
/* ------- end ---------------------------- Rayleigh.c -------------- */

