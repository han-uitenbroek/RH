/* ------- file: -------------------------- zeeman.c ----------------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Tue Jan 12 10:17:35 2010 --

       --------------------------                      ----------RH-- */

/* --- Routines to calculate Zeeman splitting patterns and relative
       strengths of Zeeman components from the labels of atomic levels.
       --                                              -------------- */

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "error.h"
#include "inputs.h"
#include "statistics.h"


/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern InputData input;
extern char messageStr[];


/* ------- begin -------------------------- determinate.c ----------- */

bool_t determinate(char *label, double g, int *n, double *S, int *L,
		   double *J)
{
  const char routineName[] = "determinate";

  char multiplet[ATOM_LABEL_WIDTH+1], *ptr, **words, orbit[3];
  int  count, multiplicity, length;

  /* --- Get the principal, spin, orbital, and angular quantum numbers
         from the atomic label --                      -------------- */

  strcpy(multiplet, label);
  ptr = multiplet + (strlen(multiplet) - 1);
  while ((*ptr != 'E')  &&  (*ptr != 'O')  &&  (ptr > multiplet))  ptr--;
  if (ptr > multiplet)
    *(ptr + 1) = '\0';
  else {
    sprintf(messageStr, "Cannot determine parity of atomic level %s", label);
    Error(WARNING, routineName, messageStr);
    return FALSE;
  }

  words = getWords(multiplet, " ", &count);
  sscanf(words[count-2], "%d", n);
  length = strlen(words[count-1]);
  sscanf(words[count-1] + length-3, "%d%s", &multiplicity, orbit);
  free(words);

  /* --- Spin quantum number --                        -------------- */

  *S = (multiplicity - 1) / 2.0;

  /* --- Orbital quantum number --                     -------------- */

  *L = getOrbital(orbit[0]);

  /* --- Total angular momentum --                     -------------- */

  *J = (g - 1.0) / 2.0;

  /* --- Composite level: cannot determine quantum numbers -- ------- */

  if (*J > *L + *S) return FALSE;

  return TRUE;
}
/* ------- end ---------------------------- determinate.c ----------- */

/* ------- begin -------------------------- getWords.c -------------- */

char **getWords(char *label, char *separator, int *count)
{
  char **theWords;
  int    length = strlen(label);

  /* --- Get the separate words that constitute the label -- -------- */

  *count = 1;
  theWords = (char **) malloc((length/2 + 1) * sizeof(char *));
  theWords[0] = strtok(label, separator);
  while ((theWords[*count] = strtok(NULL, separator)))
    *count += 1;

  return (char **) realloc(theWords, *count * sizeof(char *));
}
/* ------- end ---------------------------- getWords.c -------------- */

/* ------- begin -------------------------- effectiveLande.c -------- */

double effectiveLande(AtomicLine *line)
{
  bool_t result = TRUE;
  int    L_l, L_u, n;
  double S_l, S_u, J_l, J_u, g_l, g_u;
  Atom  *atom = line->atom;

  /* --- Determine effective Landee factor.

    See: J. Stenflo 1994, in "Solar Magnetic Fields", p. 110
         --                                            -------------- */

  if (line->g_Lande_eff != 0.0) return line->g_Lande_eff;

  result &= determinate(atom->label[line->i], atom->g[line->i],
			&n, &S_l, &L_l, &J_l);
  result &= determinate(atom->label[line->j], atom->g[line->j],
			&n, &S_u, &L_u, &J_u);

  if (result) {
    g_l = Lande(S_l, L_l, J_l);
    g_u = Lande(S_u, L_u, J_u);

    return 0.5*(g_u + g_l) + 0.25*(g_u - g_l) * 
      (J_u*(J_u + 1.0) - J_l*(J_l + 1.0));
  } else
    return 0.0;
}
/* ------- end ---------------------------- effectiveLande.c -------- */

/* ------- begin -------------------------- Lande.c ----------------- */

double Lande(double S, int L, double J)
{
  if (J == 0.0)
    return 0.0;
  else
    return 1.5 + (S*(S + 1.0) - L*(L + 1)) / (2.0*J*(J + 1.0));
}
/* ------- end ---------------------------- Lande.c ----------------- */

/* ------- begin -------------------------- ZeemanStrength.c -------- */

double ZeemanStrength(double Ju, double Mu, double Jl, double Ml)
{
  const char routineName[] = "ZeemanStrength";

  int    q, dJ;
  double s;

  /* --- Return the strength of Zeeman component (Ju, Mu) -> (Jl, Ml),
         where J and M are the total angular momentum and magnetic
         quantum numbers of the upper and lower level of a Zeeman split
         transition. These strengths are not normalised and need to be
         explicitly normalised to one using the rule:

           \Sum_{Mu, Ml} S(Ju, Mu, Jl, Ml) = 1

    See: J. Stenflo 1994, in "Solar Magnetic Fields", p. 108
         --                                            -------------- */

  q  = (int) (Ml - Mu);
  dJ = (int) (Ju - Jl);

  switch (dJ) {
  case 0:
    switch (q) {
    case  0: s = 2.0 * SQ(Mu);  break;
    case -1: s = (Ju + Mu) * (Ju - Mu + 1.0);  break;
    case  1: s = (Ju - Mu) * (Ju + Mu + 1.0);  break;
    }
    break;

  case 1:
    switch (q) {
    case  0: s = 2.0 * (SQ(Ju) - SQ(Mu));  break;
    case -1: s = (Ju + Mu) * (Ju + Mu - 1.0);  break;
    case  1: s = (Ju - Mu) * (Ju - Mu - 1.0);  break;
    }
    break;

  case -1:
    switch (q) {
    case  0: s = 2.0 * (SQ(Ju + 1.0) - SQ(Mu));  break;
    case -1: s = (Ju - Mu + 1.0) * (Ju - Mu + 2.0);  break;
    case  1: s = (Ju + Mu + 1.0) * (Ju + Mu + 2.0);  break;
    }
    break;

  default:
    sprintf(messageStr, "Invalid dJ: %d", dJ);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }   
  return s;
}
/* ------- end ---------------------------- ZeemanStrength.c -------- */

/* ------- begin -------------------------- Zeeman.c ---------------- */

ZeemanMultiplet* Zeeman(AtomicLine *line)
{
  register int n;

  bool_t result = TRUE;
  int    Ll, Lu, nq;
  double Sl, Su, Jl, Ju, Mu, Ml, norm[3], gLu, gLl;
  Atom  *atom = line->atom;
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

  if (line->g_Lande_eff != 0.0) {

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
      zm->shift[n] = zm->q[n] * line->g_Lande_eff;
    }
  } else {

    /* --- Anomalous Zeeman splitting. First get the quantum numbers
           from the labels of lower and upper level -- -------------- */

    result &= determinate(atom->label[line->i], atom->g[line->i],
			  &nq, &Sl, &Ll, &Jl);
    result &= determinate(atom->label[line->j], atom->g[line->j],
			  &nq, &Su, &Lu, &Ju);

    /* --- Count the number of components --           -------------- */

    zm->Ncomponent = 0;
    for (Ml = -Jl;  Ml <= Jl;  Ml++) {
      for (Mu = -Ju;  Mu <= Ju;  Mu++)
	if (fabs(Mu - Ml) <= 1.0) zm->Ncomponent++;
    }
    zm->q        = (int *) malloc(zm->Ncomponent * sizeof(int));
    zm->strength = (double *) malloc(zm->Ncomponent * sizeof(double));
    zm->shift    = (double *) malloc(zm->Ncomponent * sizeof(double));

    for (n = 0;  n < 3;  n++) norm[n] = 0.0;

    /* --- Fill the structure and normalize the strengths -- -------- */

    gLl = Lande(Sl, Ll, Jl);
    gLu = Lande(Su, Lu, Ju);

    n = 0;
    for (Ml = -Jl;  Ml <= Jl;  Ml++) {
      for (Mu = -Ju;  Mu <= Ju;  Mu++) {
	if (fabs(Mu - Ml) <= 1.0) {
	  zm->q[n]        = (int) (Ml - Mu);
	  zm->shift[n]    = gLl*Ml - gLu*Mu;
          zm->strength[n] = ZeemanStrength(Ju, Mu, Jl, Ml);
	  
	  norm[zm->q[n]+1] += zm->strength[n];
          n++;
	}
      }
    }
    for (n = 0;  n < zm->Ncomponent;  n++)
      zm->strength[n] /= norm[zm->q[n]+1];
  }

  return zm;
}
/* ------- end ---------------------------- Zeeman.c ---------------- */

/* ------- begin -------------------------- adjustStokesMode.c ------ */

void adjustStokesMode()
{
  register int kr, nact;

  enum StokesMode oldMode = input.StokesMode;
  Atom  *atom;
  AtomicLine *line;

  /* --- Reset Stokes mode so that the full Stokes equations are
         solved in case of input.StokesMode == FIELD_FREE, or
         input.StokesMode == POLARIZATION_FREE. --    -------------- */

  if (atmos.Nactiveatom == 0 || atmos.Stokes == FALSE ||
      input.StokesMode == NO_STOKES ||
      input.StokesMode == FULL_STOKES) return;
  else
    input.StokesMode = FULL_STOKES;

  if (oldMode == FIELD_FREE) {
    getCPU(2, TIME_START, NULL);

    for (nact = 0;  nact < atmos.Nactiveatom;  nact++) {
      atom = atmos.activeatoms[nact];

      /* --- Recalculate the profiles of polarized lines in this case */

      for (kr = 0;  kr < atom->Nline;  kr++) {
	line = &atom->line[kr];
	if (line->polarizable) {

	  /* --- First free up the space used in field-free
                 calculation --                        -------------- */
	  
	  if (!input.limit_memory) freeMatrix((void **) line->phi);
	  free(line->wphi);

	  Profile(line);
	}
      }
    }
    getCPU(2, TIME_POLL, "adjustStokesMode");
  }
}
/* ------- end ---------------------------- adjustStokesMode.c ------ */

/* ------- begin -------------------------- freeZeeman.c ------------ */

void freeZeeman(ZeemanMultiplet *zm)
{
  if (zm->Ncomponent > 0) {
    free(zm->q);
    free(zm->shift);
    free(zm->strength);
  }
}
/* ------- end ---------------------------- freeZeeman.c ------------ */

/* ------- begin -------------------------- getOrbital.c ------------ */

int getOrbital(char orbit)
{
  const char routineName[] = "getOrbital";

  int L;

  switch (orbit) {
  case 'S': L = 0;  break;
  case 'P': L = 1;  break;
  case 'D': L = 2;  break;
  case 'F': L = 3;  break;
  case 'G': L = 4;  break;
  case 'H': L = 5;  break;
  case 'I': L = 6;  break;
  case 'J': L = 7;  break;
  case 'K': L = 8;  break;
  case 'L': L = 9;  break;
  case 'M': L = 10; break;
  case 'N': L = 11;  break;
  case 'O': L = 12;  break;
  case 'Q': L = 13;  break;
  case 'R': L = 14;  break;
  case 'T': L = 15;  break;
  case 'U': L = 16;  break;
  case 'V': L = 17;  break;
  case 'W': L = 18;  break;
  case 'X': L = 19;  break;
  default: 
    sprintf(messageStr, "Invalid orbital: %c", orbit);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }
  return L;
}
/* ------- end ---------------------------- getOrbital.c ------------ */
