/* ------- file: -------------------------- readmolecule.c ----------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Fri May  7 16:57:07 2021 --

       --------------------------                      ----------RH-- */

/* --- Reads molecular data file.

       In case of an active molecule the populations of the vibrational
       levels are determined in Non-LTE, while the rotational levels
       within each vibrational level are assumed to be given by their
       respective Boltzmann factors --                 -------------- */

#include <ctype.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "spectrum.h"
#include "background.h"
#include "constant.h"
#include "error.h"
#include "inputs.h"
#include "statistics.h"


#define COMMENT_CHAR   "#"
#define N_MAX_CONFIG   10

/* --- Special case: C^13 fraction for CO --           -------------- */

#define C_13    0.011


/* --- Special case: TI fractions for TiO --           -------------- */

#define TI_46   0.08
#define TI_47   0.0743
#define TI_48   0.738
#define TI_49   0.055
#define TI_50   0.054


/* --- Function prototypes --                          -------------- */

int countTokens(char *line, char *separator);
int mrt_ascend(const void *v1, const void *v2);
int stringcompare(const void *s1, const void *s2);
char *getMoleculeID(char *molecule_file);


/* --- Global variables --                             -------------- */

extern InputData input;
extern Atmosphere atmos;
extern char messageStr[];


/* ------- begin -------------------------- readMolecule.c ---------- */

void readMolecule(Molecule *molecule, char *fileName, bool_t active)
{
  const char routineName[] = "readMolecule";
  register int k, m, n, kr, la, i;

  char    inputLine[MAX_LINE_SIZE], fitStr[20], line_data[MAX_LINE_SIZE],
          elementID[ATOM_ID_WIDTH+1], *token,
          popsFile[MAX_LINE_SIZE], labelStr[MAX_LINE_SIZE];
  bool_t  exit_on_EOF, match;
  int     Nread, Nrequired, checkPoint, J_max, v_max, Ji, Jj, nconfig;
  double  vtherm, dlambda, dlambda_wing, lambda;
  FILE   *fp_molecule;
  MolecularLine *mrt;

  getCPU(3, TIME_START, NULL);

  /* --- Open the data file for current molecule --    -------------- */

  if ((fp_molecule = fopen(fileName, "r")) == NULL) {
    sprintf(messageStr, "Unable to open inputfile %s", fileName);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  } else {
    sprintf(messageStr, " -- reading input file: %s %s", 
	    fileName, (active) ? "(active)\n\n" : "(passive)\n");


    Error(MESSAGE, routineName, messageStr);
    initMolecule(molecule);
  }
  /* --- Read molecule ID --                           -------------- */
 
  getLine(fp_molecule, COMMENT_CHAR, inputLine, exit_on_EOF=TRUE);
  Nread = sscanf(inputLine, "%s", molecule->ID);
  checkNread(Nread, Nrequired=1, routineName, checkPoint=1);

  getLine(fp_molecule, COMMENT_CHAR, inputLine, exit_on_EOF=TRUE);
  Nread = sscanf(inputLine, "%d", &molecule->charge);
  checkNread(Nread, Nrequired=1, routineName, checkPoint=2);

  /* --- We only allow molecules with 0 or +1 charge -- ------------- */

  if (molecule->charge < 0  ||  molecule->charge > 1) {
    sprintf(messageStr,
	    "Only neutral or singly charged positive molecules allowed"
	    " Molecule: %s", molecule->ID);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }

  /* --- Read list of constituent atoms --            -------------- */

  getLine(fp_molecule, COMMENT_CHAR, inputLine, exit_on_EOF=TRUE);
  molecule->Nelement = countTokens(inputLine, " ,");
  if (molecule->Nelement == 0) {
    sprintf(messageStr, "Could not identify constituents for molecule %s\n",
	    molecule->ID);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }
  molecule->pt_index =
    (int *) malloc(molecule->Nelement * sizeof(int));
  molecule->pt_count =
    (int *) malloc(molecule->Nelement * sizeof(int));
  molecule->Nnuclei = 0;

  token = strtok(inputLine, " ,");
  for (n = 0;  n < molecule->Nelement;  n++) {
    Nread = sscanf(token, "%d%s", molecule->pt_count+n, elementID);
    if (Nread == 0) {
      molecule->pt_count[n] = 1;
      Nread = sscanf(token, "%s", elementID);
    }
    UpperCase(elementID);
    if (strlen(elementID) == 1)  strcat(elementID, " ");

    match = FALSE;
    for (m = 0;  m < atmos.Nelem;  m++) {
      if (strstr(elementID, atmos.elements[m].ID)) {
	if (atmos.elements[m].abundance_set) {
	  molecule->pt_index[n] = m;
	  match = TRUE;
	}
	break;
      }
    }
    if (!match) {
      sprintf(messageStr, "Could not find element %s of molecule %s"
              " in periodic table, or element abundance not specified\n",
	      elementID, molecule->ID);
      Error(ERROR_LEVEL_2, routineName, messageStr);
    }
    molecule->Nnuclei += molecule->pt_count[n];
    token = strtok(NULL, " ,");
  }

  /* --- Read dissociation energy in eV, convert to Joule -- ------- */

  getLine(fp_molecule, COMMENT_CHAR, inputLine, exit_on_EOF=TRUE);
  Nread = sscanf(inputLine, "%lf", &molecule->Ediss);
  checkNread(Nread, Nrequired=1, routineName, checkPoint=3);
  molecule->Ediss *= EV;

  /* --- Read identification string for type of temperature fit for
         partition function and molecular equilibrium constant -- -- */

  getLine(fp_molecule, COMMENT_CHAR, inputLine, exit_on_EOF=TRUE);
  Nread = sscanf(inputLine, "%s", fitStr);
  checkNread(Nread, Nrequired=1, routineName, checkPoint=4);
  if (strstr(fitStr, "KURUCZ_70"))
    molecule->fit = KURUCZ_70;
  else if (strstr(fitStr, "KURUCZ_85"))
    molecule->fit = KURUCZ_85;
  else if (strstr(fitStr, "SAUVAL_TATUM_84"))
    molecule->fit = SAUVAL_TATUM_84;
  else if (strstr(fitStr, "IRWIN_81"))
    molecule->fit = IRWIN_81;
  else if (strstr(fitStr, "TSUJI_73"))
    molecule->fit = TSUJI_73;
  else {
    sprintf(messageStr, "Unable to regognize keyword %s in file %s",
	    fitStr, fileName);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }
  /* --- Charged molecules can only be in KURUCZ_{70,85} format -- -- */

  if (molecule->charge > 0) {
    if (molecule->fit != KURUCZ_70  &&  molecule->fit != KURUCZ_85) {
      sprintf(messageStr,
	      "Positively charged molecules MUST be in format:\n"
	      " KURUCZ_70 or\n KURUCZ_85. Molecule: %s", molecule->ID);
      Error(ERROR_LEVEL_2, routineName, messageStr);
    }
  }
  /* --- Read minimum and maximum formation temperature -- ---------- */

  getLine(fp_molecule, COMMENT_CHAR, inputLine, exit_on_EOF=TRUE);
  Nread = sscanf(inputLine, "%lf %lf", &molecule->Tmin, &molecule->Tmax);
  checkNread(Nread, Nrequired=2, routineName, checkPoint=5);

  /* --- Get fit parameters for partition function --  -------------- */

  getLine(fp_molecule, COMMENT_CHAR, inputLine, exit_on_EOF=TRUE);
  Nread = sscanf(strtok(inputLine, " "), "%d", &molecule->Npf);
  if (molecule->Npf > 0) {
    molecule->pf_coef = (double *) malloc(molecule->Npf * sizeof(double));
    for (n = molecule->Npf-1;  n >= 0;  n--)
      Nread += sscanf(strtok(NULL, " "), "%lf", molecule->pf_coef+n);
    checkNread(Nread, molecule->Npf+1, routineName, checkPoint=6);
  }
  molecule->pf = (double *) malloc(atmos.Nspace * sizeof(double));

  /* --- Get fit parameters for equilibrium constant -- ------------- */

  getLine(fp_molecule, COMMENT_CHAR, inputLine, exit_on_EOF=TRUE);
  Nread = sscanf(strtok(inputLine, " "), "%d", &molecule->Neqc);
  if (molecule->Neqc > 0) {
    molecule->eqc_coef = (double *) malloc(molecule->Neqc * sizeof(double));
    for (n = molecule->Neqc-1;  n >= 0;  n--)
      Nread += sscanf(strtok(NULL, " "), "%lf", molecule->eqc_coef+n);
    checkNread(Nread, molecule->Neqc+1, routineName, checkPoint=7);
  }

  /* --- Get molecular weight from composing atoms and calculate
         broadening velocity --                        -------------- */

  molecule->weight = 0.0;
  for (n = 0;  n < molecule->Nelement;  n++)
    molecule->weight += molecule->pt_count[n] *
      atmos.elements[molecule->pt_index[n]].weight;

  /* --- Allocate memory for total population numbers, broadening
         velocity and partition function --            -------------- */

  molecule->n = (double *) malloc(atmos.Nspace * sizeof(double));
  molecule->vbroad = (double *) malloc(atmos.Nspace * sizeof(double));

  vtherm = 2.0*KBOLTZMANN / (AMU * molecule->weight);
  for (k = 0;  k < atmos.Nspace;  k++)
    molecule->vbroad[k] = sqrt(vtherm*atmos.T[k] + SQ(atmos.vturb[k]));

  /* --- Collect the line transitions from separate data files -- --- */

  while (getLine(fp_molecule, COMMENT_CHAR, inputLine,
		 exit_on_EOF=FALSE) != EOF) {
    Nread = sscanf(inputLine, "%s", line_data);
    readMolecularLines(molecule, line_data);
  }
  /* --- Sort transitions within each molecule in wavelength -- ----- */

  qsort(molecule->mrt, molecule->Nrt, sizeof(MolecularLine), mrt_ascend);

  if (active) {
    if (molecule->Nrt == 0) {
      sprintf(messageStr, " No lines specified for active molecule %s\n",
	      molecule->ID);
      Error(ERROR_LEVEL_2, routineName, messageStr);
    }
    molecule->active = TRUE;

    /* --- This particular implementation is for the case where the
           rotational levels in each vibrational level are assumed to
           have LTE ratio's. --                        -------------- */

    v_max   = 0;
    J_max   = 0;
    nconfig = 0;

    molecule->configs    = (char *) malloc(N_MAX_CONFIG * sizeof(char));
    molecule->configs[0] = molecule->mrt[0].configi[0];
    molecule->mrt[0].ecnoi = nconfig++;

    for (kr = 0;  kr < molecule->Nrt;  kr++) {
      mrt = molecule->mrt + kr;

      v_max = MAX(v_max, MAX(mrt->vi, mrt->vj));

      Ji = (int) ((mrt->gi - 1) / 2);
      Jj = (int) ((mrt->gj - 1) / 2);
      J_max = MAX(J_max, MAX(Ji, Jj));

      /* --- Collect the unique electronic configurations --      -- */

      match = FALSE;
      for (i = 0;  i < nconfig;  i++) {
	if (mrt->configi[0] == molecule->configs[i]) {
	  match = TRUE;
          break;
	}
      }
      if (!match) {
        if (nconfig == N_MAX_CONFIG) {
	  sprintf(messageStr, "nconfig = %d > N_MAX_CONFIG", nconfig + 1);
	  Error(ERROR_LEVEL_2, routineName, messageStr);
	}
	molecule->configs[nconfig] = mrt->configi[0];
        mrt->ecnoi = nconfig++;
      }
      match = FALSE;
      for (i = 0;  i < nconfig;  i++) {
	if (mrt->configj[0] == molecule->configs[i]) {
	  match = TRUE;
          break;
	}
      }
      if (!match) {
        if (nconfig == N_MAX_CONFIG) {
	  sprintf(messageStr, "nconfig = %d > N_MAX_CONFIG", nconfig + 1);
	  Error(ERROR_LEVEL_2, routineName, messageStr);
	}
	molecule->configs[nconfig] = mrt->configj[0];
        mrt->ecnoj = nconfig++;
      }
    }
    molecule->NJ = J_max + 1;
    molecule->Nv = v_max + 1;
    molecule->Nconfig = nconfig;

    molecule->configs = (char *) realloc(molecule->configs,
				  (molecule->Nconfig + 1) * sizeof(char));
    molecule->configs[molecule->Nconfig] = '\0';
    sprintf(messageStr,
	    " --- Found %d electronic configurations for molecule %s: %s\n\n",
	    molecule->Nconfig, molecule->ID, molecule->configs);
    Error(MESSAGE, routineName, messageStr);    

    /* --- Allocate memory for the populations of each vibrational
           state in each electronic configuration. Within each of these
           the rotational states will be assumed to be in LTE -- --- */

    molecule->nv     = matrix_double(molecule->Nv * molecule->Nconfig,
				     atmos.Nspace);
    molecule->nvstar = matrix_double(molecule->Nv * molecule->Nconfig,
				     atmos.Nspace);
    molecule->pfv    = matrix_double(molecule->Nv * molecule->Nconfig,
				     atmos.Nspace);

    sprintf(popsFile, "pops_mol.%s.out", molecule->ID);
    molecule->popsFile =
      (char *) malloc((strlen(popsFile) + 1) * sizeof(char));
    strcpy(molecule->popsFile, popsFile);

    /* --- In this case vibrational partition functions for the molecule
           and each of its vibrational states are calculated in routine
           LTEmolecule --                              -------------- */ 

    LTEmolecule(molecule);

  } else {
    if (molecule->Npf > 0) {

      /* --- Otherwise, the total partition function is interpolated
	     from a polynomial in T --                 -------------- */

      for (k = 0;  k < atmos.Nspace;  k++)
	molecule->pf[k] = partfunction(molecule, atmos.T[k]);
    }
  }
  /* --- Get wavelength grid for each line if active -- ------------- */

  if (molecule->active) {
    for (kr = 0;  kr < molecule->Nrt;  kr++) {
      mrt = &molecule->mrt[kr];
      dlambda_wing = mrt->qwing * atmos.vmicro_char / CLIGHT;
      dlambda = 2.0*dlambda_wing * mrt->lambda0 / (mrt->Nlambda - 1);

      mrt->lambda = (double *) malloc(mrt->Nlambda * sizeof(double));

      lambda = mrt->lambda0 * (1.0 - dlambda_wing);
      for (la = 0;  la < mrt->Nlambda;  la++) {
	mrt->lambda[la] = lambda;
	lambda += dlambda;
      }
    }
    /* --- Allocate space for thread dependent quantities -- -------- */

    molecule->rhth =
      (rhthread *) malloc(input.Nthreads * sizeof(rhthread));
  }

  fclose(fp_molecule);
  sprintf(labelStr, "Read %s %3s",
	  (molecule->active) ? "Active" : "Molecule", molecule->ID);
  getCPU(3, TIME_POLL, labelStr);
}
/* ------- end ---------------------------- readMolecule.c ---------- */

/* ------- begin -------------------------- initMolecule.c ---------- */

void initMolecule(Molecule *molecule)
{
  molecule->popsFile = NULL;
  molecule->active = FALSE;
  molecule->fit = KURUCZ_70;
  molecule->pt_index = molecule->pt_count = 0;
  molecule->Nelement = molecule->Nnuclei = 0;
  molecule->Npf = molecule->Neqc = molecule->Nrt = 0;
  molecule->charge = 0;
  molecule->Nv = molecule->NJ = 0;
  molecule->Nconfig = 0;
  molecule->Ediss = molecule->Tmin = molecule->Tmax = molecule->weight = 0.0;
  molecule->vbroad = NULL;
  molecule->pf_coef = molecule->eqc_coef = molecule->pf = NULL;
  molecule->pfv = NULL;
  molecule->initial_solution = UNKNOWN;
  molecule->n = NULL;
  molecule->nv = molecule->nvstar = NULL;
  molecule->C_ul = NULL;
  molecule->Gamma = NULL;
  molecule->mrt = NULL;
  molecule->Ng_nv = NULL;
}
/* ------- end ---------------------------- initMolecule.c ---------- */

/* ------- begin -------------------------- freeMolecule.c ---------- */

void freeMolecule(Molecule *molecule)
{
  register int kr;

  free(molecule->pt_index);
  free(molecule->pt_count);

  if (molecule->Npf > 0) free(molecule->pf_coef);
  free(molecule->pf);
  if (molecule->Neqc > 0) free(molecule->eqc_coef);

  free(molecule->n);
  free(molecule->vbroad);

  if (molecule->active) {
    free(molecule->configs);
    freeMatrix((void **) molecule->nv);
    freeMatrix((void **) molecule->nvstar);
    freeMatrix((void **) molecule->pfv);

    free(molecule->popsFile);
    free(molecule->rhth);
  }

  if (molecule->Nrt > 0) {
    for (kr = 0;  kr < molecule->Nrt;  kr++)
      freeMolLine(&molecule->mrt[kr]);
    free(molecule->mrt);
  }
}
/* ------- end ---------------------------- freeMolecule.c ---------- */

/* ------- begin -------------------------- initMolLine.c ----------- */

void initMolLine(MolecularLine *mrt, enum type line_type)
{
  /* --- Initialize transition data structure. Set pointers to NULL - */

  mrt->type = line_type;
  mrt->configi[0] = mrt->configj[0] = '\0';
  mrt->symmetric = TRUE;
  mrt->Voigt = TRUE;
  mrt->vi = mrt->vj = 0; 
  mrt->Nlambda = mrt->Nblue = 0;
  mrt->lambda0 = mrt->isotope_frac = 0.0;
  mrt->lambda = NULL;
  mrt->phi = NULL;
  mrt->wphi = NULL;
  mrt->Ei = mrt->Ej = mrt->gi = mrt->gj = 0.0;
  mrt->Grad = mrt->qwing = mrt->qcore = 0.0;
  mrt->Aji = mrt->Bji = mrt->Bij = 0.0;

  mrt->polarizable = FALSE;
  mrt->molecule = NULL;
  mrt->zm = NULL;
  mrt->g_Lande_eff = 0.0;
}
/* ------- end ---------------------------- initMolLine.c ----------- */

/* ------- begin -------------------------- freeMolLine.c ----------- */

void freeMolLine(MolecularLine *mrt)
{
  /* --- Free allocated memory for active transition structure line - */

  Molecule *molecule = mrt->molecule;
  
  if (molecule->active) {
    freeMatrix((void **) mrt->phi);
    free(mrt->phi);
  }
  
  if (atmos.Stokes &&
      input.StokesMode == FULL_STOKES &&
      mrt->polarizable) {
    freeZeeman(mrt->zm);
  }
}
/* ------- end ---------------------------- freeMolLine.c ----------- */

/* ------- begin -------------------------- mrt_ascend.c ------------ */
 
int mrt_ascend(const void *v1, const void *v2)
{
  double f1, f2;

  /* --- Used for sorting molecular transitions by wavelength -- ---- */

  f1 = ((MolecularLine *) v1)->lambda0;
  f2 = ((MolecularLine *) v2)->lambda0;
  if (f1 < f2)
    return -1;
  else if (f1 > f2)
    return 1;
  else
    return 0;
}
/* ------- end ---------------------------- mrt_ascend.c ------------ */

/* ------- begin -------------------------- stringcompare.c --------- */

int stringcompare(const void *s1, const void *s2)
{
  /* --- Wrapper around strcmp to get the proper prototype for the
         comparing function (ie., with void rather than char
         parameters). --                               -------------- */

  return strcmp((const char *) s1, (const char *) s2);
}
/* ------- end ---------------------------- stringcompare.c --------- */

/* ------- begin -------------------------- readMolecularLines.c ---- */

/* --- Reads molecular line data and computes molecular opacities.
 
       Recognized formats:

         -- GOORVITCH94:

          CO line data in the format of D. Goorvitch 1994, ApJS 95, 535
          (vibration-rotation transitions of the X configuration).


         -- KURUCZ_CD18:

          Molecular line data from Bob Kurucz's CD 18

  203.6264 -7.917  2.5    83.924  2.5 -49177.701 108X00F1   A07E1   16
|    |    |    |    |    |    |    |    |    |    |    |    |    |    |    |
0    5   10        20        30        40        50        60        70

   The numbers mean:

     203.6264  --  Wavelength [nm]
      -7.917   --  log(gf) value
       2.5     --  J first level (usually lower level)
      83.924   --  Energy first level (- indicates predicted, + observed)
                   in units of [cm^-1]
       2.5     --  J second level (usually upper level)
  -49177.701   --  Energy second level [cm^-1]
    108X00F1   --  Molecule ID (OH in this case) plus ID first level:
                   X configuration, vibration level 0, F parity 
       A07E1   --  ID second level: A configuration, vibration level 7,
                   parity E
        16     --  Isotope ID (16 refers to O16 in this case)


         -- KURUCZ_TIO

          Molecular line list for TiO molecules

  708.9002 -1.213 68.0  9920.304 67.0 24022.776822a04p1   b08m1   46
|    |    |    |    |    |    |    |    |    |    |    |    |    |    |    |
0    5   10        20        30        40        50        60        70

   Meaning of the numbers as above.

   Parity translation p --> E, and m --> F, for the orbital angular
   momentum pointing parallel or anti-parallel to the internuclear axis.

       --                                              -------------- */

void readMolecularLines(struct Molecule *molecule, char *line_data)
{
  const char routineName[] = "readMolecularLines";
  register int kr;

  char   inputLine[MAX_LINE_SIZE], branchStr[2], type_string[25],
         format_string[15], Lambdai[2], Lambdaj[2], Hundi[2], Hundj[2];
  bool_t exit_on_EOF;
  enum   type line_type;
  int    Nrt, Nread, Nrequired, checkPoint, isotope, Nlambda, Ji, Jj;
  double dipole, C, strength, gf, qwing, lambda_air, log_gf, lambda0;
  FILE  *fp_lines;
  MolecularLine *mrt;

  C = 2.0*PI * (Q_ELECTRON/EPSILON_0) * (Q_ELECTRON/M_ELECTRON) / CLIGHT;

  /* --- Open the data file --                         -------------- */
 
  if ((fp_lines = fopen(line_data, "r")) == NULL) {
    sprintf(messageStr, "Unable to open inputfile %s", line_data);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  } else {
    sprintf(messageStr, " -- reading input file: %s\n", line_data);
    Error(MESSAGE,routineName, messageStr);
  }
  /* --- Get count, type, and format of the lines in this file -- -- */

  getLine(fp_lines, COMMENT_CHAR, inputLine, exit_on_EOF=TRUE);
  Nread = sscanf(inputLine, "%d %s %s", &Nrt,
		 type_string, format_string);
  checkNread(Nread, Nrequired=3, routineName, checkPoint=1);

  if (strstr(type_string, "VIBRATION_ROTATION"))
    line_type = VIBRATION_ROTATION;
  else if (strstr(type_string, "MOLECULAR_ELECTRONIC"))
    line_type = MOLECULAR_ELECTRONIC;
  else {
    sprintf(messageStr, "Unknown molecular line type: %s", type_string);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }
  molecule->mrt = (MolecularLine *)
    realloc(molecule->mrt, (molecule->Nrt + Nrt) * sizeof(MolecularLine));

  getLine(fp_lines, COMMENT_CHAR, inputLine, exit_on_EOF=TRUE);
  Nread = sscanf(inputLine, "%d %lf", &Nlambda, &qwing);
  checkNread(Nread, Nrequired=2, routineName, checkPoint=2);

  /* --- Read line data line by line --                -------------- */

  if (strstr(format_string, "GOORVITCH94")) {
    for (kr = molecule->Nrt;  kr < molecule->Nrt + Nrt;  kr++) {
      mrt = molecule->mrt + kr;
      initMolLine(mrt, VIBRATION_ROTATION);

      mrt->molecule = molecule;
      
      strcpy(mrt->configi, "X");
      strcpy(mrt->configj, "X");

      mrt->Nlambda = Nlambda;
      mrt->qwing   = qwing;

      getLine(fp_lines, COMMENT_CHAR, inputLine, exit_on_EOF=TRUE);
      Nread = sscanf(inputLine, "%lf %lf %lf %lf %lf %lf %d %d %s %d %d",
		     &mrt->lambda0, &dipole, &mrt->Aji, &mrt->Ei, &gf,
		     &strength, &mrt->vj, &mrt->vi, branchStr, &Ji,
		     &isotope);
      checkNread(Nread, Nrequired=11, routineName, checkPoint=3);
    
      mrt->lambda0 = (CM_TO_M / NM_TO_M) / mrt->lambda0;
      mrt->Ei *= (HPLANCK * CLIGHT) / CM_TO_M;
      mrt->Ej  = mrt->Ei + HPLANCK*CLIGHT / (mrt->lambda0 * NM_TO_M);

      mrt->gi = 2*Ji + 1;
      switch (branchStr[0]) {
      case 'P':
	Jj = Ji - 1;
	break;
      case 'Q':
	Jj = Ji;
	break;
      case 'R':
	Jj = Ji + 1;
	break;
      default:
	sprintf(messageStr, "Invalid value of branch string in file %s: %s\n"
		" Line: %f [nm]\n  Valid options are P, R and Q",
		line_data, branchStr, mrt->lambda0);
	Error(ERROR_LEVEL_2, routineName, messageStr);
      }
      mrt->gj = 2*Jj + 1;

      mrt->Bji = CUBE(mrt->lambda0*NM_TO_M) /
	(2.0*HPLANCK*CLIGHT) * mrt->Aji;
      mrt->Bij = (mrt->gj / mrt->gi) * mrt->Bji;
    
      /* --- Temporary measure to accommodate $^13$CO --   ---------- */  

      switch (isotope) {
      case 12:
      case 26: mrt->isotope_frac = 1.0 - C_13;  break;
      case 13:
      case 36: mrt->isotope_frac = C_13;        break;
      default: mrt->isotope_frac = 1.0;
      }
      
      mrt->polarizable = FALSE;
    }
  } else if (strstr(format_string, "KURUCZ_CD18") ||
	     strstr(format_string, "KURUCZ_NEW")  ||
	     strstr(format_string, "KURUCZ_TIO")) {

    C = 2 * PI * (Q_ELECTRON/EPSILON_0) * (Q_ELECTRON/M_ELECTRON) / CLIGHT;

    for (kr = molecule->Nrt;  kr < molecule->Nrt + Nrt;  kr++) {
      mrt = &molecule->mrt[kr];
      initMolLine(mrt, line_type);

      mrt->molecule = molecule;
      mrt->Nlambda  = Nlambda;
      mrt->qwing    = qwing;
      mrt->isotope_frac = 1.0;

      getLine(fp_lines, COMMENT_CHAR, inputLine, exit_on_EOF=TRUE);

      if (strstr(format_string, "KURUCZ_TIO")) {
	Nread = sscanf(inputLine,
		       "%10lf %7lf %5lf %10lf %5lf %10lf",
		       &lambda_air, &log_gf,
		       &mrt->gi, &mrt->Ei, &mrt->gj, &mrt->Ej);
      } else {
	Nread = sscanf(inputLine,
		       "%10lf %7lf %5lf %10lf %5lf %11lf",
		       &lambda_air, &log_gf,
		       &mrt->gi, &mrt->Ei, &mrt->gj, &mrt->Ej);
      }
      checkNread(Nread, Nrequired=6, routineName, checkPoint=3);
    
      mrt->Ei = (HPLANCK * CLIGHT) / CM_TO_M * fabs(mrt->Ei);
      mrt->Ej = (HPLANCK * CLIGHT) / CM_TO_M * fabs(mrt->Ej);
      mrt->gi = 2*mrt->gi + 1;
      mrt->gj = 2*mrt->gj + 1;

      lambda0  = (HPLANCK * CLIGHT) / (mrt->Ej - mrt->Ei);
      mrt->Aji = C / SQ(lambda0) * POW10(log_gf)/mrt->gj;
      mrt->Bji = CUBE(lambda0) / (2.0 * HPLANCK * CLIGHT) * mrt->Aji;
      mrt->Bij = (mrt->gj / mrt->gi) * mrt->Bji;
      mrt->lambda0 = lambda0 / NM_TO_M;

      if (strstr(format_string, "KURUCZ_CD18")) {

	/* --- Electronic configuration --             -------------- */

	mrt->configi[0] = inputLine[52];
	mrt->configj[0] = inputLine[60];
	mrt->configi[1] = '\0';
	mrt->configj[1] = '\0';

	/* --- Parity designation --                   -------------- */

	mrt->parityi[0] = inputLine[55];
	mrt->parityj[0] = inputLine[63];
	mrt->parityi[1] = '\0';
	mrt->parityj[1] = '\0';
  

	/* --- Vibrational quantum numbers --            ------------ */

	Nread = sscanf(inputLine+53, "%2d", &mrt->vi);
	Nread = sscanf(inputLine+61, "%2d", &mrt->vj);

	/* --- Subbranch (F1, F2, ...., E1, E2, ...etc) -- ---------- */

	Nread = sscanf(inputLine+56, "%1d", &mrt->subi);
	Nread = sscanf(inputLine+64, "%1d", &mrt->subj);

      } else if (strstr(format_string, "KURUCZ_NEW")) {

	/* --- Electronic configuration --             -------------- */

	strncpy(mrt->configi, inputLine+52, 2);
	strncpy(mrt->configj, inputLine+60, 2);
	mrt->configi[2] = '\0';
	mrt->configj[2] = '\0';

	/* --- Parity designation --                   -------------- */

	mrt->parityi[0] = inputLine[56];
	mrt->parityj[0] = inputLine[64];
	mrt->parityi[1] = '\0';
	mrt->parityj[1] = '\0';
  
	/* --- Vibrational quantum numbers --            ------------ */

	Nread = sscanf(inputLine+54, "%2d", &mrt->vi);
	Nread = sscanf(inputLine+62, "%2d", &mrt->vj);

	/* --- Subbranch (F1, F2, ...., E1, E2, ...etc) -- ---------- */

	Nread = sscanf(inputLine+57, "%1d", &mrt->subi);
	Nread = sscanf(inputLine+65, "%1d", &mrt->subj);

      } else if (strstr(format_string, "KURUCZ_TIO")) {

	/* --- Electronic configuration --             -------------- */

	strncpy(mrt->configi, inputLine+50, 2);
	strncpy(mrt->configj, inputLine+59, 2);
	mrt->configi[2] = '\0';
	mrt->configj[2] = '\0';

	/* --- Parity designation --                   -------------- */

	mrt->parityi[0] = ((inputLine[53] == 'p') ? 'E' : 'F');
	mrt->parityj[0] = ((inputLine[61] == 'p') ? 'E' : 'F');
	mrt->parityi[1] = '\0';
	mrt->parityj[1] = '\0';
  
	/* --- Vibrational quantum numbers --            ------------ */

	Nread = sscanf(inputLine+51, "%2d", &mrt->vi);
	Nread = sscanf(inputLine+59, "%2d", &mrt->vj);

	/* --- Subbranch (F1, F2, ...., E1, E2, ...etc) -- ---------- */

	Nread = sscanf(inputLine+54, "%1d", &mrt->subi);
	Nread = sscanf(inputLine+62, "%1d", &mrt->subj);

        /* --- Isotope fraction for Ti --                ------------ */

        Nread = sscanf(inputLine+66, "%2d", &isotope);
        switch (isotope) {
	case 46: mrt->isotope_frac = TI_46;  break;
	case 47: mrt->isotope_frac = TI_47;  break;
	case 48: mrt->isotope_frac = TI_48;  break;
	case 49: mrt->isotope_frac = TI_49;  break;
	case 50: mrt->isotope_frac = TI_50;  break;
	default: mrt->isotope_frac = 1.0;
	}
      }
      /* --- read additional data for Zeeman polarization -- -------- */

      if (strlen(inputLine) > 71) {

	Nread = sscanf(inputLine+71, "%1s %1s %lf %1s %1s %lf",
		       Hundi, Lambdai, &mrt->Si,
		       Hundj, Lambdaj, &mrt->Sj);

	/* --- Determine the coupling case according to Hund -- ----- */

        switch (Hundi[0]) {
	case 'A': mrt->Hundi = CASE_A;  break;
	case 'B': mrt->Hundi = CASE_B;  break;
        default:
          sprintf(messageStr, "Unsupported Hund's case: %s\n", Hundi);
	  Error(ERROR_LEVEL_2, routineName, messageStr);
	}
        switch (Hundj[0]) {
	case 'A': mrt->Hundj = CASE_A;  break;
	case 'B': mrt->Hundj = CASE_B;  break;
        default:
          sprintf(messageStr, "Unsupported Hund's case: %s\n", Hundj);
	  Error(ERROR_LEVEL_2, routineName, messageStr);
	}
	/* --- Orbital angular momemtum Lambda along nuclear axis

               S --> Sigma
               P --> Pi
               D --> Delta
               F --> Phi
           --                                          -------------- */

        switch (Lambdai[0]) {
	case 'S': mrt->Lambdai = 0;  break;
	case 'P': mrt->Lambdai = 1;  break;
	case 'D': mrt->Lambdai = 2;  break;
	case 'F': mrt->Lambdai = 3;  break;
        default:
          sprintf(messageStr,
		  "Unsupported orbital projection Lambda: %s\n", Lambdai);
	  Error(ERROR_LEVEL_2, routineName, messageStr);
	}
        switch (Lambdaj[0]) {
	case 'S': mrt->Lambdaj = 0;  break;
	case 'P': mrt->Lambdaj = 1;  break;
	case 'D': mrt->Lambdaj = 2;  break;
	case 'F': mrt->Lambdaj = 3;  break;
        default:
          sprintf(messageStr,
		  "Unsupported orbital projection Lambda: %s\n", Lambdaj);
	  Error(ERROR_LEVEL_2, routineName, messageStr);
	}
	/* --- Set the polarized flag only if the atmosphere has magnetic
	       fields --                               -------------- */
	
	if (atmos.Stokes) mrt->polarizable = TRUE;
      } else
	mrt->polarizable = FALSE;
    }
  } else {
    sprintf(messageStr, "Unknown molecular line format: %s", format_string);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }
  molecule->Nrt += Nrt;

  fclose(fp_lines);
  sprintf(messageStr, " --- read %d %s lines for molecule %2s\n\n",
	  Nrt, type_string, molecule->ID);
  Error(MESSAGE, routineName, messageStr);
}
/* ------- end ---------------------------- readMolecularLines.c ---- */

/* ------- begin -------------------------- readMolecularModels.c --- */

void readMolecularModels(void)
{
  const char routineName[] = "readMolecularModels";
  register int n, m, i;

  char    filename[MAX_LINE_SIZE], inputLine[MAX_LINE_SIZE],
          actionKey[MAX_KEYWORD_SIZE], popsKey[MAX_KEYWORD_SIZE],
          popsFile[MAX_LINE_SIZE], *moleculeID;
  bool_t  active, exit_on_EOF;
  int     Nread, Nrequired, checkPoint;
  FILE   *fp_molecules;
  Molecule *molecule;
  Element *element;

  getCPU(2, TIME_START, NULL);

  /* --- Open input file for molecular models --       -------------- */

  if ((fp_molecules = fopen(input.molecules_input, "r")) == NULL) {
    sprintf(messageStr, "Unable to open input file %s",
	    input.molecules_input);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }
  /* --- Get the number of molecular models to be read -- ----------- */

  getLine(fp_molecules, COMMENT_CHAR, inputLine, exit_on_EOF=TRUE);
  Nread = sscanf(inputLine, "%d", &atmos.Nmolecule);
  checkNread(Nread, Nrequired=1, routineName, checkPoint=1);

  atmos.molecules =
    (Molecule *) malloc(atmos.Nmolecule * sizeof(Molecule));

  atmos.H2 = NULL;
  atmos.OH = NULL;
  atmos.CH = NULL;

  /* --- Read molecular data for the various molecules -- ----------- */

  for (n = 0;  n < atmos.Nmolecule;  n++) {

    getLine(fp_molecules, COMMENT_CHAR, inputLine, exit_on_EOF=TRUE);
    Nread = sscanf(inputLine, "%s %s %s %s",
		   filename, actionKey, popsKey, popsFile);
    checkNread(Nread, Nrequired=3, routineName, checkPoint=2);

    moleculeID = getMoleculeID(filename);

    /* --- The first molecule MUST be H2 --            -------------- */

    if (n == 0  &&  !strstr(moleculeID, "H2")) {
      sprintf(messageStr,
	      "\n First molecule should be H2\n");
      Error(ERROR_LEVEL_2, routineName, messageStr);
    }
    /* --- Check for duplicate molecule IDs --         -------------- */

    for (m = 0;  m < n;  m++) {
      if (strcmp(moleculeID, atmos.molecules[m].ID) == 0) {
	sprintf(messageStr, "Duplicate model for molecule %s\n",
		moleculeID);
	Error(ERROR_LEVEL_2, routineName, messageStr);
      }
    }

    molecule = &atmos.molecules[n];
    readMolecule(molecule, filename,
		 active=(strstr(actionKey, "ACTIVE") ? TRUE : FALSE));

    /* --- Create an alias to the hydrogen molecule  structure -- ----- */

    if (strcmp(molecule->ID, "H2") == 0) atmos.H2 = molecule;

    /* --- Store pointers to OH and CH molecules, if appropriate,
           for use in bf opacity calculation (see ohchbf.c) -- ------ */                      

    if (strcmp(molecule->ID, "OH") == 0) atmos.OH = molecule;
    if (strcmp(molecule->ID, "CH") == 0) atmos.CH = molecule;
  
    /* --- Set flag for initial soltion --             -------------- */

    if (strstr(popsKey, "OLD_POPULATIONS")) {
      molecule->initial_solution = OLD_POPULATIONS;
      if (Nread < 4) {
	sprintf(messageStr,
		"No file with OLD_POPULATIONS specified for molecule: %s\n",
		moleculeID);
	Error(ERROR_LEVEL_2, routineName, messageStr);
      }
      molecule->popsFile =
	(char *) malloc((strlen(popsFile) + 1) * sizeof(char));
      strcpy(molecule->popsFile, popsFile);
    } else if (strstr(popsKey, "LTE_POPULATIONS")) {
      molecule->initial_solution = LTE_POPULATIONS;
    }

    /* --- If popsKey is not recognized --             -------------- */

    if (molecule->initial_solution == UNKNOWN) {
      sprintf(messageStr,
	      "Unknown initial solution specified for molecule: %s\n"
              "Has to be LTE_POPULATIONS for molecules\n",
	      moleculeID);
      Error(ERROR_LEVEL_2, routineName, messageStr);
    } 

  }
  fclose(fp_molecules);
  
  /* --- Figure out for each element in which molecule it may be
         bound and store the indices of those molecules -- ---------- */

  for (n = 0;  n < atmos.Nelem;  n++) {
    element = &atmos.elements[n];
    element->Nmolecule = 0;
    element->mol_index = (int *) malloc(atmos.Nmolecule * sizeof(int));

    for (m = 0;  m < atmos.Nmolecule;  m++) {
      molecule = &atmos.molecules[m];
      for (i = 0;  i < molecule->Nelement;  i++) {
	if (molecule->pt_index[i] == n) {
	  element->mol_index[element->Nmolecule++] = m;
	  break;
	}
      }
    }
    if (element->Nmolecule == 0) {
      free(element->mol_index);
      element->mol_index = NULL;
    } else
      element->mol_index = (int *)
	realloc(element->mol_index, element->Nmolecule * sizeof(int));
  }
  /* --- Allocate memory for H^- --                    -------------- */

  atmos.nHmin = (double *) malloc(atmos.Nspace * sizeof(double));

  getCPU(2, TIME_POLL, "Read molecular input");
}
/* ------- end ---------------------------- readMolecularModels.c --- */

/* ------- begin -------------------------- getMoleculeID.c --------- */

char *getMoleculeID(char *molecule_file)
{
  const char routineName[] = "getMoleculeID";
  register int n;
  static char moleculeID[MOLECULE_ID_WIDTH + 1];

  char   inputLine[MAX_LINE_SIZE];
  bool_t exit_on_EOF;
  FILE  *fp_molecule;

  if ((fp_molecule = fopen(molecule_file, "r")) == NULL) {
    sprintf(messageStr, "Unable to open inputfile %s", molecule_file);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }

  /* --- Read molecule ID --                           -------------- */
 
  getLine(fp_molecule, COMMENT_CHAR, inputLine, exit_on_EOF=TRUE);
  sscanf(inputLine, "%s", moleculeID);
  for (n = 0;  n < (int) strlen(moleculeID);  n++)
    moleculeID[n] = toupper(moleculeID[n]);

  fclose(fp_molecule);
  return moleculeID;
}
/* ------- end ---------------------------- getMoleculeID.c --------- */


/* ------- begin -------------------------- countTokens.c ----------- */

int countTokens(char *line, char *separator)
{
  char *string, *token;
  int   count = 0;

  /* --- Count the number of tokens separated by separator in line -- */

  string = (char *) malloc((strlen(line) + 1) * sizeof(char));
  strcpy(string, line);

  token = strtok(string, separator);
  while (token) {
    count++;
    token = strtok(NULL, separator);
  }

  free(string);
  return count;
}
/* ------- end ---------------------------- countTokens.c ----------- */
