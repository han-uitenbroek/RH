/* ------- file: -------------------------- barklem.c ---------------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Fri May 14 19:18:30 2021 --

       --------------------------                      ----------RH-- */

/* --- Routines to deal with Collisional broadening as formulated by
       Barklem et al.

  See: - Anstee & O'Mara 1995, MNRAS 276, 859-866 (s-p, p-s)
       - Barklem & O'Mara 1997, MNRAS 290, 102 (p-d, d-p)
       - Barklem, O'Mara & Ross 1998, MNRAS 296, 1057-1060 (d-f, f-d)
       - Barklem, O'Mara 1998, MNRAS 300, 863-871
       --                                              -------------- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "constant.h"
#include "error.h"


#define BARKLEM_SP_DATA     "../../Atoms/Barklem_spdata.dat"
#define BARKLEM_SP_NS       21
#define BARKLEM_SP_NP       18
#define BARKLEM_SP_NEFF1    1.0
#define BARKLEM_SP_NEFF2    1.3

#define BARKLEM_PD_DATA     "../../Atoms/Barklem_pddata.dat"
#define BARKLEM_PD_NP       18
#define BARKLEM_PD_ND       18
#define BARKLEM_PD_NEFF1    1.3
#define BARKLEM_PD_NEFF2    2.3

#define BARKLEM_DF_DATA     "../../Atoms/Barklem_dfdata.dat"
#define BARKLEM_DF_ND       18
#define BARKLEM_DF_NF       18
#define BARKLEM_DF_NEFF1    2.3
#define BARKLEM_DF_NEFF2    3.3

#define BARKLEM_DELTA_NEFF  0.1


/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern char messageStr[];


/* ------- begin -------------------------- readBarklemTable.c ------ */

bool_t readBarklemTable(enum Barklemtype type, Barklemstruct *bs)
{
  register int n, i, j;
  const char routineName[] = "readBarklemTable";

  char    filename[MAX_LINE_SIZE], inputLine[MAX_LINE_SIZE], *charptr;
  int     nread;
  double  neff1_0, neff2_0;
  FILE   *fp_Barklem;

  switch (type) {
  case SP:
    strcpy(filename, BARKLEM_SP_DATA);
    bs->N1 = BARKLEM_SP_NS;
    bs->N2 = BARKLEM_SP_NP;

    neff1_0 = BARKLEM_SP_NEFF1;
    neff2_0 = BARKLEM_SP_NEFF2;
    break;

  case PD:
    strcpy(filename, BARKLEM_PD_DATA);
    bs->N1 = BARKLEM_PD_NP;
    bs->N2 = BARKLEM_PD_ND;

    neff1_0 = BARKLEM_PD_NEFF1;
    neff2_0 = BARKLEM_PD_NEFF2;
    break;

  case DF:
    strcpy(filename, BARKLEM_DF_DATA);
    bs->N1 = BARKLEM_DF_ND;
    bs->N2 = BARKLEM_DF_NF;

    neff1_0 = BARKLEM_DF_NEFF1;
    neff2_0 = BARKLEM_DF_NEFF2;
    break;
  }

  if ((fp_Barklem = fopen(filename, "r")) == NULL) {
    sprintf(messageStr, "Unable to open input file %s", filename);
    Error(ERROR_LEVEL_1, routineName, messageStr);
    return FALSE;
  }

  bs->neff1 = (double *) malloc(bs->N1 * sizeof(double));
  for (n = 0;  n < bs->N1;  n++)
    bs->neff1[n] = neff1_0 + n * BARKLEM_DELTA_NEFF;
 
  bs->neff2 = (double *) malloc(bs->N2 * sizeof(double));
  for (n = 0;  n < bs->N2;  n++)
    bs->neff2[n] = neff2_0 + n * BARKLEM_DELTA_NEFF;

  bs->cross = matrix_double(bs->N1, bs->N2);
  bs->alpha = matrix_double(bs->N1, bs->N2);

  for (n = 0;  n < 3;  n++)
    charptr = fgets(inputLine, MAX_LINE_SIZE, fp_Barklem);

  for (i = 0;  i < bs->N1;  i++)
    for (j = 0;  j < bs->N2;  j++) {
      nread = fscanf(fp_Barklem, "%lf", &bs->cross[i][j]);
  }
  for (n = 0;  n < 2;  n++)
    charptr = fgets(inputLine, MAX_LINE_SIZE, fp_Barklem);

  for (i = 0;  i < bs->N1;  i++)
    for (j = 0;  j < bs->N2;  j++) {
      nread = fscanf(fp_Barklem, "%lf", &bs->alpha[i][j]);
  }

  fclose(fp_Barklem);
  return TRUE;
}
/* ------- end ---------------------------- readBarklemTable.c------- */

/* ------- begin -------------------------- getBarklemcross.c ------- */

bool_t getBarklemcross(Barklemstruct *bs, RLK_Line *rlk)
{
  const char routineName[] = "getBarklemcross";

  int index;
  double Z, neff1, neff2, findex1, findex2, reducedmass, meanvelocity,
         crossmean, E_Rydberg, deltaEi, deltaEj;
  Element *element;

  element = &atmos.elements[rlk->pt_index - 1];

  /* --- Note: ABO tabulations are valid only for neutral atoms -- -- */

  if (rlk->stage > 0)
    return FALSE;

  if ((deltaEi = element->ionpot[rlk->stage] - rlk->level_i.E) <= 0.0)
    return FALSE;
  if ((deltaEj = element->ionpot[rlk->stage] - rlk->level_j.E) <= 0.0)
    return FALSE;

  Z = (double) (rlk->stage + 1);
  E_Rydberg = E_RYDBERG / (1.0 + M_ELECTRON / (element->weight * AMU));
  neff1 = Z * sqrt(E_Rydberg / deltaEi);
  neff2 = Z * sqrt(E_Rydberg / deltaEj);

  if (rlk->level_i.L > rlk->level_j.L) SWAPDOUBLE(neff1, neff2);
  
  if (neff1 < bs->neff1[0] || neff1 > bs->neff1[bs->N1-1])
    return FALSE;
  Locate(bs->N1, bs->neff1, neff1, &index);
  findex1 =
    (double) index + (neff1 - bs->neff1[index]) / BARKLEM_DELTA_NEFF;

  if (neff2 < bs->neff2[0] || neff2 > bs->neff2[bs->N2-1])
    return FALSE;
  Locate(bs->N2, bs->neff2, neff2, &index);
  findex2 =
    (double) index + (neff2 - bs->neff2[index]) / BARKLEM_DELTA_NEFF;

  /* --- Find interpolation in table --                -------------- */

  rlk->cross = cubeconvol(bs->N2, bs->N1,
			  bs->cross[0], findex2, findex1);
  rlk->alpha = cubeconvol(bs->N2, bs->N1,
			  bs->alpha[0], findex2, findex1);

  reducedmass  = AMU / (1.0/atmos.H->weight + 1.0/element->weight);
  meanvelocity = sqrt(8.0 * KBOLTZMANN / (PI * reducedmass));
  crossmean    = SQ(RBOHR) * pow(meanvelocity / 1.0E4, -rlk->alpha);

  rlk->cross *= 2.0 * pow(4.0/PI, rlk->alpha/2.0) * 
    exp(gammln((4.0 - rlk->alpha)/2.0)) * meanvelocity * crossmean;  

  rlk->vdwaals = BARKLEM;
  return TRUE;
}
/* ------- end ---------------------------- getBarklemcross.c ------- */

/* ------- begin -------------------------- getBarklemactivecross.c - */

bool_t getBarklemactivecross(AtomicLine *line)
{
  bool_t determined = TRUE, useBarklem = FALSE;
  int index, Ll, Lu, nq, i, j, ic;
  double Sl, Su, Jl, Ju;
  double Z, neff1, neff2, findex1, findex2, reducedmass, meanvelocity,
         crossmean, E_Rydberg, deltaEi, deltaEj;
  Atom *atom;
  Barklemstruct bs;
 
  atom = line->atom;
  j = line->j;
  i = line->i;

  /* --- ABO tabulations are only valid for neutral atoms -- -------- */

  if (atom->stage[i] > 0)
    return FALSE;

  /* --- Get the quantum numbers for orbital angular momentum -- ---- */

  determined &= determinate(atom->label[i], atom->g[i],
			    &nq, &Sl, &Ll, &Jl);
  determined &= determinate(atom->label[j], atom->g[j],
			    &nq, &Su, &Lu, &Ju);

  /* --- See if one of the Barklem cases applies --    -------------- */

  if (determined) {
    if ((Ll == S_ORBIT && Lu == P_ORBIT) ||
	(Ll == P_ORBIT && Lu == S_ORBIT)) {
      useBarklem = readBarklemTable(SP, &bs);
    } else if ((Ll == P_ORBIT && Lu == D_ORBIT) ||
	       (Ll == D_ORBIT && Lu == P_ORBIT)) {
      useBarklem = readBarklemTable(PD, &bs);
    } else if ((Ll == D_ORBIT && Lu == F_ORBIT) ||
	       (Ll == F_ORBIT && Lu == D_ORBIT)) {
      useBarklem = readBarklemTable(DF, &bs);
    }
  }
  if (!determined || !useBarklem) return FALSE;

  /* --- Determine the index of the appropriate continuum level -- -- */

  Z = atom->stage[j] + 1;
  for (ic = j + 1;  atom->stage[ic] < atom->stage[j]+1;  ic++);

  deltaEi   = atom->E[ic] - atom->E[i];
  deltaEj   = atom->E[ic] - atom->E[j];
  E_Rydberg = E_RYDBERG / (1.0 + M_ELECTRON / (atom->weight * AMU));

  neff1 = Z * sqrt(E_Rydberg / deltaEi);
  neff2 = Z * sqrt(E_Rydberg / deltaEj);

  if (Ll > Lu) SWAPDOUBLE(neff1, neff2);
  
  /* --- Interpolate according to effective principal quantum number  */

  if (neff1 < bs.neff1[0] || neff1 > bs.neff1[bs.N1-1])
    return FALSE;
  Locate(bs.N1, bs.neff1, neff1, &index);
  findex1 =
    (double) index + (neff1 - bs.neff1[index]) / BARKLEM_DELTA_NEFF;

  if (neff2 < bs.neff2[0] || neff2 > bs.neff2[bs.N2-1])
    return FALSE;
  Locate(bs.N2, bs.neff2, neff2, &index);
  findex2 =
    (double) index + (neff2 - bs.neff2[index]) / BARKLEM_DELTA_NEFF;

  /* --- Find interpolation in table --                -------------- */

  line->cvdWaals[0] = cubeconvol(bs.N2, bs.N1,
				 bs.cross[0], findex2, findex1);
  line->cvdWaals[1] = cubeconvol(bs.N2, bs.N1,
				 bs.alpha[0], findex2, findex1);

  reducedmass  = AMU / (1.0/atmos.atoms[0].weight + 1.0/atom->weight);
  meanvelocity = sqrt(8.0 * KBOLTZMANN / (PI * reducedmass));
  crossmean    = SQ(RBOHR) * pow(meanvelocity / 1.0E4, -line->cvdWaals[1]);

  line->cvdWaals[0] *= 2.0 * pow(4.0/PI, line->cvdWaals[1]/2.0) * 
    exp(gammln((4.0 - line->cvdWaals[1])/2.0)) * meanvelocity * crossmean;  

  /* --- Use UNSOLD for the contribution of Helium atoms -- ---------- */

  line->cvdWaals[2] = 1.0;
  line->cvdWaals[3] = 0.0;

  return TRUE;
}
/* ------- end ---------------------------- getBarklemactivecross.c -- */
