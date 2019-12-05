/* ------- file: -------------------------- broad.c -----------------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Fri Nov  9 11:07:55 2007 --

       --------------------------                      ----------RH-- */

/* --- Routines for line broadening, including van der Waals and
       Stark broadening by electronic collisions. --   -------------- */
 
#include <ctype.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "constant.h"
#include "error.h"

/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern char messageStr[];


/* ------- begin -------------------------- VanderWaals.c ----------- */

void VanderWaals(AtomicLine *line, double *GvdW)
{
  /* --- Compute van der Waals broadening in Lindholm theory with
         Unsold's approximation for the interaction coefficient C_6.

    See: Traving 1960, "Uber die Theorie der Druckverbreiterung
           von Spektrallinien", p 91-97
      -- Mihalas 1978, p. 282ff, and Table 9-1.

         Gamma = 8.08 * vrel^3/5 * C_6^2/5 * atmos.H.ntotal_neutral


         Or with parametrization using Smirnov-Roueff potential.

    See: DeRidder & van Rensbergen 1976, A&A Suppl. 23, 147-165

         Gamma = alpha * T^{beta} * atmos.H.ntotal_neutral 


         Or with the Anstee, Barklem & O'Mara formalism.

    See: Anstee & O'Mara 1995, MNRAS 276, 859-866
         Barklem & O'Mara 1998, MNRAS 300, 863-871

         Gamma =
          (4/pi)^(alpha/2) Gam((4-alpha)/2) v_0 sigma(v_0)(vmean/v_0)^(1-alpha)

                                                           -- ------- */

  const char routineName[] = "VanderWaals";
  register int k;

  int     i, j, ic, Z;
  double  vrel35_H, vrel35_He, fourPIeps0, deltaR, cross, gammaCorrH,
          gammaCorrHe, *T = atmos.T, C625;
  Atom    *atom = line->atom;
  Element *He = &atmos.elements[1];

  j = line->j;
  i = line->i;

  if (line->vdWaals == UNSOLD  || line->vdWaals == BARKLEM) {
    fourPIeps0 = 4.0 * PI * EPSILON_0;

    vrel35_He = pow(8.0*KBOLTZMANN/(PI * AMU * atom->weight) * 
		    (1.0 + atom->weight/He->weight), 0.3);

    Z = atom->stage[j] + 1;
    for (ic = j + 1;  atom->stage[ic] < atom->stage[j]+1;  ic++);

    deltaR = SQ(E_RYDBERG/(atom->E[ic] - atom->E[j])) -
      SQ(E_RYDBERG/(atom->E[ic] - atom->E[i]));
    C625  = pow(2.5 * (SQ(Q_ELECTRON)/fourPIeps0) * (ABARH/fourPIeps0) *
		2*PI * SQ(Z*RBOHR)/HPLANCK * deltaR, 0.4);
  }

  switch (line->vdWaals) {
  case UNSOLD:

    /* --- Relative velocity of radiator and perturber with Maxwellian
           velocity distributions --                   -------------- */

    vrel35_H  = pow(8.0*KBOLTZMANN/(PI * AMU * atom->weight) * 
		    (1.0 + atom->weight/atmos.H->weight), 0.3);

    cross = 8.08 * (line->cvdWaals[0]*vrel35_H +
		    line->cvdWaals[2]*He->abund*vrel35_He) * C625;

    for (k = 0;  k < atmos.Nspace;  k++)
      GvdW[k] = cross * pow(T[k], 0.3);

    break;

  case RIDDER_RENSBERGEN:

    /* --- alpha = 1.0E-8 * cvdW[0]  (Hydrogen broadening)
                 = 1.0E-9 * cvdW[1]  (Helium broadening) -- --------- */

    gammaCorrH  = 1.0E-8 * CUBE(CM_TO_M) *
      pow(1.0 + atmos.H->weight/atom->weight, line->cvdWaals[1]);
    gammaCorrHe = 1.0E-9 * CUBE(CM_TO_M) *
      pow(1.0 + He->weight/atom->weight, line->cvdWaals[3]);

    for (k = 0;  k < atmos.Nspace;  k++)
      GvdW[k] = gammaCorrH*line->cvdWaals[0] *
	pow(T[k], line->cvdWaals[1]) + 
	gammaCorrHe*line->cvdWaals[2] *
	pow(T[k], line->cvdWaals[3]) * He->abund;

    break;

  case BARKLEM:

    /* --- Use UNSOLD for the Helium contribution --   -------------- */

    cross = 8.08 * line->cvdWaals[2]*He->abund*vrel35_He * C625;

    for (k = 0;  k < atmos.Nspace;  k++)
      GvdW[k] = line->cvdWaals[0] *
	pow(atmos.T[k], (1.0 - line->cvdWaals[1])/2.0) +
	cross * pow(T[k], 0.3);

    break;

  default:
    sprintf(messageStr,
	    "Unknown method for van der Waals broadening %d, line %d -> %d",
	    line->vdWaals, line->j, line->i);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }
  /* --- Multiply with the Hydrogen ground level population -- ------ */

  for (k = 0;  k < atmos.Nspace;  k++) GvdW[k] *= atmos.H->n[0][k];
}
/* ------- end ---------------------------- VanderWaals.c ----------- */

/* ------- begin -------------------------- Stark.c ----------------- */

#define AVERAGE_ATOMIC_WEIGHT  28.0

void Stark(AtomicLine *line, double *GStark)
{
  /* --- Quadratic Stark broadening by electrons and singly charged 
         ions. 

         Gamma = 11.37 * vrel^1/3 * C_4^2/3 * (ne + nion)

	 Use estimate for C_4 from Traving.
     
    See: Traving 1960, "Uber die Theorie der Druckverbreiterung
                 von Spektrallinien", p 93

      -- Mihalas 1978, p. 282ff, and Table 9-1.

      -- David F. Gray, Observation and Analysis of Stellar
         Photospheres (1992), 2nd ed., p. 216, eq. 11.33


         if line->cStark < 0 then Gamma = abs(line->cStark) * ne
         --                                            -------------- */

  const char routineName[] = "Stark";
  register int k, ic;

  int    Z;
  double C4, C, Cm, m_electron = M_ELECTRON/AMU, cStark23, cStark,
         neff_u, neff_l, m_atom_avg = AVERAGE_ATOMIC_WEIGHT, vrel,
         E_Rydberg;
  Atom  *atom = line->atom;

  if (line->cStark < 0.0) {
    cStark = fabs(line->cStark);
    for (k = 0;  k < atmos.Nspace;  k++)
      GStark[k] = cStark * atmos.ne[k];
  } else {
    /* --- Constants for relative velocity. We assume that nion = ne
           (see Gray), and that the average atomic weight of ionic
           perturbers is given by AVERAGE_ATOMIC_WEIGHT. -- --------- */
    
    C  = 8.0 * KBOLTZMANN / (PI * AMU * atom->weight);
    Cm = pow(1.0 + atom->weight/m_electron, 0.16666667) + 
      pow(1.0 + atom->weight/m_atom_avg, 0.16666667);
    
    /* --- Find core charge Z and effective quantum numbers neff_u
           and neff_l for upper and lower level --     -------------- */
    
    Z = atom->stage[line->i] + 1;
    for (ic = line->i + 1;
	 ((atom->stage[ic] < atom->stage[line->i]+1) &&
	  (ic < atom->Nlevel));  ic++);
    if (atom->stage[ic] == atom->stage[line->i]) {
      sprintf(messageStr, "Cannot find overlying continuum for level %d",
	      line->i);
      Error(ERROR_LEVEL_2, routineName, messageStr);
    }
    E_Rydberg = E_RYDBERG / (1.0 + M_ELECTRON / (atom->weight * AMU));
    neff_l = Z * sqrt(E_Rydberg / (atom->E[ic] - atom->E[line->i]));
    neff_u = Z * sqrt(E_Rydberg / (atom->E[ic] - atom->E[line->j]));
    
    C4 = (SQ(Q_ELECTRON) / (4.0 * PI * EPSILON_0)) * RBOHR *
      (2.0*PI * SQ(RBOHR) / HPLANCK) / (18.0 * SQ(Z)*SQ(Z)) *
	(SQ(neff_u*(5.0*SQ(neff_u) + 1.0)) -
	 SQ(neff_l*(5.0*SQ(neff_l) + 1.0)));
    cStark23 = 11.37 * pow(line->cStark * C4, 0.66666667);
    
    for (k = 0;  k < atmos.Nspace;  k++) {
      vrel = pow(C * atmos.T[k], 0.16666667) * Cm;
      GStark[k] = cStark23 * vrel * atmos.ne[k];
    }
  }
}
/* ------- end ---------------------------- Stark.c ----------------- */

/* ------- begin -------------------------- StarkLinear.c ----------- */

void StarkLinear(AtomicLine *line, double *GStark)
{
  /* --- Linear Stark broadening by electrons for hydrogen lines.
     
    See: K. Sutton (1978), JQSRT 20, 333-343

         GStark = a_1 * [0.60 * (n_u^2 - n_l^2) * (N_e)^(2/3) * CM_TO_M^2]
         --                                            -------------- */

  const char routineName[] = "StarkLinear";
  register int k;

  char   config[4], *ptr;
  int    n_upper, n_lower;
  double a1, C;
  Atom  *atom = line->atom;

  if (strstr(atom->ID, "H ") == NULL) {
    sprintf(messageStr, "Model is not a hydrogen atom: %s", atom->ID);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }
  /* --- Find principal quantum number of lower and upper level -- -- */

  sscanf(atom->label[line->i], "H I %s", config);
  ptr = config;  while (isdigit(*ptr)) ptr++;  *ptr = ' ';
  sscanf(config, "%d", &n_lower);

  sscanf(atom->label[line->j], "H I %s", config);
  ptr = config;  while (isdigit(*ptr)) ptr++;  *ptr = ' ';
  sscanf(config, "%d", &n_upper);

  if (n_upper - n_lower == 1)
    a1 = 0.642;
  else
    a1 = 1.0;

  C = a1 * 0.6 * (SQ(n_upper) - SQ(n_lower)) * SQ(CM_TO_M);
  for (k = 0;  k < atmos.Nspace;  k++)
    GStark[k] = C * pow(atmos.ne[k], 0.66666667);
}
/* ------- end ---------------------------- StarkLinear.c ----------- */

/* ------- begin -------------------------- Damping.c --------------- */

void Damping(AtomicLine *line, double *adamp)
{
  register int k;

  double cDop, *Qelast;
  Atom *atom;

  cDop = (NM_TO_M * line->lambda0) / (4.0 * PI);

  atom = line->atom;
  Qelast = (double *) calloc(atmos.Nspace, sizeof(double));

  /* --- Add van der Waals broadening --             -------------- */

  if ((line->cvdWaals[0] > 0.0) || (line->cvdWaals[2] > 0.0)) {
    VanderWaals(line, adamp);
    for (k = 0;  k < atmos.Nspace;  k++) Qelast[k] += adamp[k];
  }
  /* --- Add Quadratic Stark broadening --           -------------- */

  if (line->cStark != 0.0) {
    Stark(line, adamp);
    for (k = 0;  k < atmos.Nspace;  k++) Qelast[k] += adamp[k];
  }
  /* --- Add Linear Stark broadening for hydrogen only -- --------- */

  if (strstr(atom->ID, "H ")) {
    StarkLinear(line, adamp);
    for (k = 0;  k < atmos.Nspace;  k++) Qelast[k] += adamp[k];
  }
  /* --- Store the total rate of elastic collisions in case of PRD  */

  if (line->PRD  &&  line->Qelast != NULL) {
    for (k = 0;  k < atmos.Nspace;  k++)
      line->Qelast[k] = Qelast[k];
  }
  for (k = 0;  k < atmos.Nspace;  k++)
    adamp[k] = (line->Grad + Qelast[k]) * cDop / atom->vbroad[k];

  free(Qelast);
}
/* ------- end ---------------------------- Damping.c --------------- */

/* ------- begin -------------------------- MolecularDamping.c ------ */

void MolecularDamping(MolecularLine *mrt, double *adamp)
{
  register int k;

  double cDop;
  Molecule *molecule;

  cDop = (NM_TO_M * mrt->lambda0) / (4.0 * PI);
  molecule = mrt->molecule;

  /* --- For now only natural broadening due to the line itself. -- - */

  for (k = 0;  k < atmos.Nspace;  k++)
    adamp[k] = mrt->Aji * cDop / molecule->vbroad[k];
}
/* ------- end ---------------------------- MolecularDamping.c ------ */

