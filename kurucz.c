/* ------- file: -------------------------- kurucz.c ----------------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Fri May 28 16:02:49 2021 --

       --------------------------                      ----------RH-- */

/* --- Routines to deal with a Kurucz-format line list. These are lines
       to be used as LTE background.

       The input file should contain a list of Kurucz linelist files,
       one line each, in the format described on the kurucz web page:
         http://kurucz.harvard.edu/linelists.html

FORMAT(F11.4,F7.3,F6.2,F12.3,F5.2,1X,A10,F12.3,F5.2,1X,A10,
3F6.2,A4,2I2,I3,F6.3,I3,F6.3,2I5,1X,A1,A1,1X,A1,A1,i1,A3.2I5,I6) 

 1 wavelength (nm)  air above 200 nm   F11.4
 2 log gf  F7.3
 3 element code = element number + charge/100.  F6.2
 4 first energy level in cm-1   F12.3
 5 J for first level   F5.1
   blank for legibility   1X
 6 label field for first level   A10
 7 second energy level in cm-1   F12.3
        (negative energies are predicted or extrapolated}  
 8 J for second level   F5.1
   blank for legibility   1X
 9 label field for second level   A10
10 log of radiative damping constant, Gamma Rad  F6.2 or F6.3
11 log of stark damping constant/electron number. Gamma Stark  F6.2 or F6.3
12 log of van der Waals damping constant/neutral hydrogen number, 
       Gamma van der Waals   F6.2 or F6.3
13 reference that can be expanded in subdirectory LINES   A4  
14 non-LTE level index for first level   I2
15 non-LTE level index for second level   I2
16 isotope number   I3
17 hyperfine component log fractional strength  F6.3
18 isotope number  (for diatomics there are two and no hyperfine)   I3
19 log isotopic abundance fraction   F6.3
20 hyperfine shift for first level in mK to be added to E  I5
21 hyperfine shift for second level in mK to be added to E'  I5
   the symbol "F" for legibilty   1X
22 hyperfine F for the first level    I1
23 note on character of hyperfine data for first level: z none, ? guessed  A1
   the symbol "-" for legibility    1X
24 hyperfine F' for the second level  I1
25 note on character of hyperfine data for second level: z none, ? guessed  A1
26 1-digit code, sometimes for line strength classes   I1
27 3-character code such as AUT for autoionizing    A3  
28 lande g for first level times 1000   I5
29 lande g for second level times 1000   I5
30 isotope shift of wavelength in mA 


 Note: The periodic table index of Kurucz starts at 1 (for Hydrogen).
       In the RH routines counting starts at 0 in the atmos.elements
       structure. This discrepancy is accounted for in the opacity
       calculation (rlk_opacity.c).
       --                                              -------------- */

#include <ctype.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "background.h"
#include "spectrum.h"
#include "constant.h"
#include "inputs.h"
#include "error.h"

#define COMMENT_CHAR             "#"
#define RLK_RECORD_LENGTH        160
#define Q_WING                   20.0
#define MILLI                    1.0E-03
#define ANGSTROM_TO_NM           0.1
#define MAX_GAUSS_DOPPLER        7.0
#define USE_TABULATED_WAVELENGTH 1

#define RLK_LABEL_LENGTH  10


/* --- Function prototypes --                          -------------- */

double RLKProfile(RLK_Line *rlk, int k, int mu, bool_t to_obs,
		  double lambda,
		  double *phi_Q, double *phi_U, double *phi_V,
		  double *psi_Q, double *psi_U, double *psi_V);
void   RLKZeeman(RLK_Line *rlk);
double RLKLande(RLK_level* level);
void   initRLK(RLK_Line *rlk);
bool_t RLKdet_level(char* label, RLK_level *level);
double getJK_K(char c);
void   getUnsoldcross(RLK_Line *rlk);
void   free_BS(Barklemstruct *bs);


/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern InputData input;
extern char messageStr[];


/* ------- begin -------------------------- readKuruczLines.c ------- */

void readKuruczLines(char *inputFile)
{
  const char routineName[] = "readKuruczLines";
  const double  C = 2.0*PI * (Q_ELECTRON/EPSILON_0) * 
                             (Q_ELECTRON/M_ELECTRON) / CLIGHT;

  char   inputLine[RLK_RECORD_LENGTH+1], listName[MAX_LINE_SIZE],
         filename[MAX_LINE_SIZE], Gvalues[20+1], elem_code[7],
         labeli[RLK_LABEL_LENGTH+1], labelj[RLK_LABEL_LENGTH+1],
        *commentChar = COMMENT_CHAR;
  bool_t swap_levels, determined, useBarklem;
  int    Nline, Nread, Nrequired, checkPoint, hfs_i, hfs_j, gL_i, gL_j,
    iso_dl, Li, Lj, tmp;
  double lambda0, Ji, Jj, Grad, GStark, GvdWaals, pti,
         Ei, Ej, gf, lambda_air;
  RLK_Line *rlk;
  Barklemstruct bs_SP, bs_PD, bs_DF;
  FILE  *fp_Kurucz, *fp_linelist;

  if (!strcmp(inputFile, "none")) return;

  /* --- Read in the data files for Barklem collisional broadening -- */

  readBarklemTable(SP, &bs_SP);
  readBarklemTable(PD, &bs_PD);
  readBarklemTable(DF, &bs_DF);

  labeli[RLK_LABEL_LENGTH] = '\0';
  labelj[RLK_LABEL_LENGTH] = '\0';

  if ((fp_Kurucz = fopen(inputFile, "r")) == NULL) {
    sprintf(messageStr, "Unable to open input file %s", inputFile);
    Error(ERROR_LEVEL_1, routineName, messageStr);
    return;
  }
  /* --- Go through each of the linelist files listed in input file - */  

  while (getLine(fp_Kurucz, commentChar, listName, FALSE) != EOF) {
    Nread = sscanf(listName, "%s", filename);
    if ((fp_linelist = fopen(filename, "r")) == NULL) {
      sprintf(messageStr, "Unable to open input file %s", filename);
      Error(ERROR_LEVEL_1, routineName, messageStr);
    }
    /* --- Count the number of lines in this file --   -------------- */

    Nline = 0;
    while (fgets(inputLine, RLK_RECORD_LENGTH+1, fp_linelist) != NULL)
      if (*inputLine != *commentChar) Nline++;
    rewind(fp_linelist);

    if (atmos.Nrlk == 0) atmos.rlk_lines = NULL;
    atmos.rlk_lines = (RLK_Line *)
      realloc(atmos.rlk_lines, (Nline + atmos.Nrlk) * sizeof(RLK_Line));

    /* --- Read lines from file --                     -------------- */

    rlk = atmos.rlk_lines + atmos.Nrlk;
    while (fgets(inputLine, RLK_RECORD_LENGTH+1, fp_linelist) != NULL) {
      if (*inputLine != *commentChar) {

        initRLK(rlk);

	Nread = sscanf(inputLine, "%lf %lf %s %lf",
		       &lambda_air, &gf, (char *) &elem_code, &Ei);

        /* --- Ionization stage and periodic table index -- --------- */

        sscanf(elem_code, "%d.%d", &rlk->pt_index, &rlk->stage);

	Nread += sscanf(inputLine+53, "%lf", &Ej);

	Ei = fabs(Ei) * (HPLANCK * CLIGHT) / CM_TO_M;
	Ej = fabs(Ej) * (HPLANCK * CLIGHT) / CM_TO_M;

	/* --- Beware: the Kurucz linelist has upper and lower levels
	       of a transition in random order. Therefore, we have to
               check for the lowest energy of the two and use that as
               lower level --                          -------------- */

	if (Ej < Ei) {
	  swap_levels = TRUE; 
	  rlk->level_i.E = Ej;
	  rlk->level_j.E = Ei;
	  strncpy(labeli, inputLine+69, RLK_LABEL_LENGTH);
	  strncpy(labelj, inputLine+41, RLK_LABEL_LENGTH);
	} else {
	  swap_levels = FALSE;
	  rlk->level_i.E = Ei;
          rlk->level_j.E = Ej;
	  strncpy(labeli, inputLine+41, RLK_LABEL_LENGTH);
	  strncpy(labelj, inputLine+69, RLK_LABEL_LENGTH);
	}

	Nread += sscanf(inputLine+35, "%lf", &Ji);
	Nread += sscanf(inputLine+63, "%lf", &Jj);
	if (swap_levels) SWAPDOUBLE(Ji, Jj);
	
	rlk->level_i.J = Ji;
	rlk->level_i.g = 2*Ji + 1;
	rlk->level_j.J = Jj;
	rlk->level_j.g = 2*Jj + 1;

	if (USE_TABULATED_WAVELENGTH) {
	  
	  /* --- In this case use tabulated wavelength and adjust 
	         upper-level energy --                 -------------- */

	  air_to_vacuum(1, &lambda_air, &lambda0);
	  lambda0 *= NM_TO_M;
	  rlk->level_j.E = rlk->level_i.E +  (HPLANCK * CLIGHT) / lambda0;
	} else {
	  /* --- Else use energy levels to calculate lambda0 -- ----- */
	  
	  lambda0 = (HPLANCK * CLIGHT) / (rlk->level_j.E - rlk->level_i.E);
	}
	rlk->Aji = C / SQ(lambda0) * POW10(gf) / rlk->level_j.g;
	rlk->Bji = CUBE(lambda0) / (2.0 * HPLANCK * CLIGHT) * rlk->Aji;
	rlk->Bij = (rlk->level_j.g / rlk->level_i.g) * rlk->Bji;

        /* --- Store in nm --                          -------------- */

	rlk->lambda0 = lambda0 / NM_TO_M;

	/* --- Get quantum numbers for angular momentum and spin -- - */

        determined = (RLKdet_level(labeli, &rlk->level_i) &&
		      RLKdet_level(labelj, &rlk->level_j));
        rlk->polarizable = (atmos.Stokes && determined);

        /* --- Line broadening --                      -------------- */

	strncpy(Gvalues, inputLine+79, 20); // JdlCR: read 3 more characters from the vdW label 
	Nread += sscanf(Gvalues, "%lf %lf %lf", &Grad, &GStark, &GvdWaals);

	if (GStark != 0.0) 
	  rlk->GStark = POW10(GStark) * CUBE(CM_TO_M);
	else
	  rlk->GStark = 0.0;

	if (GvdWaals != 0.0){
	  if(GvdWaals < 20.0) // JdlCR: If larger than 20, we are providing sigma.alpha
	    rlk->GvdWaals = POW10(GvdWaals) * CUBE(CM_TO_M);
	  else{ // User provided Barklem cross sections, integer part is sigma and decimal part is alpha
	    tmp = (int)GvdWaals;
	    rlk->cross = (double)tmp;    // integer part
	    rlk->alpha = GvdWaals - tmp; // decimal part
	  }
	}else
	  rlk->GvdWaals = 0.0;

        /* --- If possible use Barklem formalism --    -------------- */
	
	useBarklem = FALSE;
	if (determined &&
	    rlk->level_i.cpl == LS_COUPLING &&
	    rlk->level_j.cpl == LS_COUPLING) {

	  /* --- JdlCR: changed to lower "l" value, the search will
	     fail if "l" is not provided as both will be zero --- */
	  
	  Li = rlk->level_i.l;
	  Lj = rlk->level_j.l;
	  
	  if ((Li == S_ORBIT && Lj == P_ORBIT) ||
              (Li == P_ORBIT && Lj == S_ORBIT)) {
	    useBarklem = getBarklemcross(&bs_SP, rlk);
	  } else if ((Li == P_ORBIT && Lj == D_ORBIT) ||
		     (Li == D_ORBIT && Lj == P_ORBIT)) {
	    useBarklem = getBarklemcross(&bs_PD, rlk);
	  } else if ((Li == D_ORBIT && Lj == F_ORBIT) ||
		     (Li == F_ORBIT && Lj == D_ORBIT)) {
	    useBarklem = getBarklemcross(&bs_DF, rlk);
	  }
	}
	/* --- Else use good old Unsoeld --            -------------- */

        if (!useBarklem) {
	  getUnsoldcross(rlk);
	}
	/* --- Radiative broadening --                 -------------- */

	if (Grad != 0.0) {
	  rlk->Grad = POW10(Grad);
	} else {

	  /* --- Just take the Einstein Aji value--    -------------- */     

	  rlk->Grad = rlk->Aji;
	}
	/* --- Isotope and hyperfine fractions and slpittings -- ---- */

	Nread += sscanf(inputLine+106, "%d", &rlk->isotope);
	Nread += sscanf(inputLine+108, "%lf", &rlk->isotope_frac);
	rlk->isotope_frac = POW10(rlk->isotope_frac);
	Nread += sscanf(inputLine+117, "%lf", &rlk->hyperfine_frac);
	rlk->hyperfine_frac = POW10(rlk->hyperfine_frac);
	Nread += sscanf(inputLine+123, "%5d%5d", &hfs_i, &hfs_j);
	rlk->level_i.hfs = ((double) hfs_i) * MILLI * KBOLTZMANN;
	rlk->level_j.hfs = ((double) hfs_j) * MILLI * KBOLTZMANN;

	/* --- Effective Lande factors --              -------------- */

	Nread += sscanf(inputLine+143, "%5d%5d", &gL_i, &gL_j);
	rlk->level_i.gL = gL_i * MILLI;
	rlk->level_j.gL = gL_j * MILLI;
	if (swap_levels) {
	  SWAPDOUBLE(rlk->level_i.hfs, rlk->level_j.hfs);
	  SWAPDOUBLE(rlk->level_i.gL, rlk->level_j.gL);
	}

	/*      Nread += sscanf(inputLine+154, "%d", &iso_dl); */
	iso_dl = 0;
	rlk->iso_dl = iso_dl * MILLI * ANGSTROM_TO_NM;

	checkNread(Nread, Nrequired=17, routineName, checkPoint=1);
	rlk++;
      }
    }
    fclose(fp_linelist);

    sprintf(messageStr, "Read %d Kurucz lines from file %s\n",
	    Nline, listName);
    Error(MESSAGE, routineName, messageStr);
    atmos.Nrlk += Nline;
  }

  fclose(fp_Kurucz);

  free_BS(&bs_SP);
  free_BS(&bs_PD);
  free_BS(&bs_DF);
}
/* ------- end ---------------------------- readKuruczLines.c ------- */

/* ------- begin -------------------------- rlk_ascend.c ------------ */
 
int rlk_ascend(const void *v1, const void *v2)
{
  double lambda1 = ((RLK_Line *) v1)->lambda0,
         lambda2 = ((RLK_Line *) v2)->lambda0;

  /* --- Used for sorting transitions by wavelength -- -------------- */

  if (lambda1 < lambda2)
    return -1;
  else if (lambda1 > lambda2)
    return 1;
  else
    return 0;
}
/* ------- end ---------------------------- rlk_ascend.c ------------ */

/* ------- begin -------------------------- rlk_locate.c ------------ */

void rlk_locate(int N, RLK_Line *lines, double lambda, int *low)
{
  int  high, index, increment;

  /* --- Locate position wavelength lambda in Kurucz line list. Assume
         that lines have been sorted in order of ascending wavelength */

  if ((*low <= 0)  ||  (*low > N-1)) {

    /* --- Input guess not useful here, go to bisection --  --------- */

    *low = 0;
    high = N;
  } else {

    /* --- Else hunt up or down to bracket value --    -------------- */ 

    increment = 1;
    if (lambda >= lines[*low].lambda0) {
      high = *low + increment;
      if (*low == N-1) return;

      /* --- Hunt up --                                -------------- */

      while (lambda >= lines[high].lambda0) {
	*low = high;
	increment += increment;
	high = *low + increment;
        if (high >= N) { high = N;  break; }
      }
    } else {
      high = *low;
      if (*low == 0) return;

      /* --- Hunt down --                              -------------- */

      while (lambda <= lines[*low].lambda0) {
	high = *low;
	increment += increment;
	*low = high - increment;
        if (*low <= 0) { *low = 0;  break; }
      }
    }
  }
  /* --- Bisection algorithm --                        -------------- */

  while (high - *low > 1) {
    index = (high + *low) >> 1;
    if (lambda >= lines[index].lambda0)
      *low = index;
    else
      high = index;
  }
}
/* ------- end ---------------------------- rlk_locate.c ------------ */

/* ------- begin -------------------------- rlk_opacity.c ----------- */

flags rlk_opacity(double lambda, int nspect, int mu, bool_t to_obs,
                  double *chi, double *eta, double *scatt, double *chip)
{
  register int k, n, kr;

  bool_t contributes, hunt;
  int    Nwhite, Nblue, Nred, NrecStokes;
  double dlamb_wing, *pf, dlamb_char, hc_la, ni_gi, nj_gj, lambda0, kT,
         Bijhc_4PI, twohnu3_c2, hc, fourPI, hc_4PI,
        *eta_Q, *eta_U, *eta_V, eta_l,
        *chi_Q, *chi_U, *chi_V, chi_l, *chip_Q, *chip_U, *chip_V,
         phi, phi_Q, phi_U, phi_V, psi_Q, psi_U, psi_V,
         epsilon, C, C2_atom, C2_ion, C3, dE, x;
  Atom *metal;
  AtomicLine *line;
  Element *element;
  RLK_Line *rlk;
  flags backgrflags;

  /* --- Calculate the LTE opacity at wavelength lambda due to atomic
         transitions stored in atmos.rlk_lines --      -------------- */

  backgrflags.hasline     = FALSE;
  backgrflags.ispolarized = FALSE;

  /* --- If wavelength outside our list return without calculation -- */

  dlamb_char = lambda * Q_WING * (atmos.vmicro_char / CLIGHT);
  if (lambda < atmos.rlk_lines[0].lambda0 - dlamb_char ||
      lambda > atmos.rlk_lines[atmos.Nrlk-1].lambda0 + dlamb_char) {
   return backgrflags;
  }

  hc     = HPLANCK * CLIGHT;
  fourPI = 4.0 * PI;
  hc_4PI = hc / fourPI;

  if (input.rlkscatter) {
    C       = 2 * PI * (Q_ELECTRON/EPSILON_0) *
                (Q_ELECTRON/M_ELECTRON) / CLIGHT;
    C2_atom = 2.15E-6;
    C2_ion  = 3.96E-6;
  }

  pf = (double *) malloc(atmos.Nspace * sizeof(double));

  /* --- locate wavelength lambda in table of lines -- -------------- */

  Nwhite = 0;
  rlk_locate(atmos.Nrlk, atmos.rlk_lines, lambda, &Nwhite);
  Nblue = Nwhite;
  while (atmos.rlk_lines[Nblue].lambda0 + dlamb_char > lambda &&
	 Nblue > 0)  Nblue--;
  Nred = Nwhite;
  while (atmos.rlk_lines[Nred].lambda0 - dlamb_char < lambda &&
	 Nred < atmos.Nrlk-1)  Nred++;

  /* --- Initialize the contribution for this wavelength and angle -- */

  if (Nred >= Nblue) {
    if (atmos.Stokes) {
      NrecStokes = 4;

      /* --- Use pointers to sub-arrays for Q, U, and V -- ---------- */

      chi_Q = chi + atmos.Nspace;
      chi_U = chi + 2*atmos.Nspace;
      chi_V = chi + 3*atmos.Nspace;

      eta_Q = eta + atmos.Nspace;
      eta_U = eta + 2*atmos.Nspace;
      eta_V = eta + 3*atmos.Nspace;

      if (input.magneto_optical) {
        chip_Q = chip;
        chip_U = chip + atmos.Nspace;
        chip_V = chip + 2*atmos.Nspace;

        for (k = 0;  k < 3*atmos.Nspace;  k++) chip[k] = 0.0;
      }
    } else
      NrecStokes = 1;

    for (k = 0;  k < NrecStokes * atmos.Nspace;  k++) {
      chi[k] = 0.0;
      eta[k] = 0.0;
    }
    if (input.rlkscatter) {
      for (k = 0;  k < atmos.Nspace;  k++) scatt[k] = 0.0;
    }
  }
  /* --- Add opacities from lines at this wavelength -- ------------- */

  for (n = Nblue;  n <= Nred;  n++) {
    rlk = &atmos.rlk_lines[n];
    if (fabs(rlk->lambda0 - lambda) <= dlamb_char) {      
      element = &atmos.elements[rlk->pt_index - 1];

      /* --- Check whether partition function is present for this
	     stage, and if abundance is set --         -------------- */

      if ((rlk->stage < element->Nstage - 1) && element->abundance_set) {
	contributes = TRUE;
	if ((metal = element->model) != NULL) {

          /* --- If an explicit atomic model is present check that we
	         do not already account for this line in this way - - */

	  for (kr = 0;  kr < metal->Nline;  kr++) {
	    line = metal->line + kr;
	    dlamb_wing = line->lambda0 * line->qwing *
	      (atmos.vmicro_char / CLIGHT);
	    if (fabs(lambda - line->lambda0) <= dlamb_wing &&
		metal->stage[line->i] == rlk->stage) {
	      contributes = FALSE;
	      break;
	    }
	  }
	}
      } else
	contributes = FALSE;

      /* --- Get opacity from line --                  -------------- */

      if (contributes) {
	hc_la      = (HPLANCK * CLIGHT) / (rlk->lambda0 * NM_TO_M);
	Bijhc_4PI  = hc_4PI * rlk->Bij * rlk->isotope_frac *
	  rlk->hyperfine_frac * rlk->level_i.g;
	twohnu3_c2 = rlk->Aji / rlk->Bji;

	if (input.rlkscatter) {
	  if (rlk->stage == 0) {
	    x  = 0.68;
	    C3 = C / (C2_atom * SQ(rlk->lambda0 * NM_TO_M));
	  } else {
	    x  = 0.0;
	    C3 = C / (C2_ion * SQ(rlk->lambda0 * NM_TO_M));
	  }

	  dE = rlk->level_j.E - rlk->level_i.E;
	}
        /* --- Set flag that line is present at this wavelength -- -- */

	backgrflags.hasline = TRUE;
	if (rlk->polarizable) {
	  backgrflags.ispolarized = TRUE;
	  if (rlk->zm == NULL) RLKZeeman(rlk);
	}

        if (element->n == NULL) {
	  element->n = matrix_double(element->Nstage, atmos.Nspace);
	  LTEpops_elem(element);
	}
        Linear(atmos.Npf, atmos.Tpf, element->pf[rlk->stage],
	       atmos.Nspace, atmos.T, pf, hunt=TRUE);

	for (k = 0;  k < atmos.Nspace;  k++) {
	  phi = RLKProfile(rlk, k, mu, to_obs, lambda,
			   &phi_Q, &phi_U, &phi_V,
			   &psi_Q, &psi_U, &psi_V);

	  if (phi){
	    kT    = 1.0 / (KBOLTZMANN * atmos.T[k]);
	    ni_gi = element->n[rlk->stage][k] *
	      exp(-rlk->level_i.E * kT - pf[k]);
            nj_gj = ni_gi * exp(-hc_la * kT);

	    chi_l = Bijhc_4PI * (ni_gi - nj_gj);
	    eta_l = Bijhc_4PI * twohnu3_c2 * nj_gj;

	    if (input.rlkscatter) {
	      epsilon = 1.0 / (1.0 + C3 * pow(atmos.T[k], 1.5) /
			       (atmos.ne[k] * 
				pow(KBOLTZMANN * atmos.T[k] / dE, 1 + x)));

              scatt[k] += (1.0 - epsilon) * chi_l * phi;
	      chi_l    *= epsilon;
              eta_l    *= epsilon;
	    }

	    chi[k] += chi_l * phi;
	    eta[k] += eta_l * phi;

	    if (rlk->zm != NULL && rlk->Grad) {
	      chi_Q[k] += chi_l * phi_Q;
	      chi_U[k] += chi_l * phi_U;
	      chi_V[k] += chi_l * phi_V;

	      eta_Q[k] += eta_l * phi_Q;
	      eta_U[k] += eta_l * phi_U;
	      eta_V[k] += eta_l * phi_V;

	      if (input.magneto_optical) {
		chip_Q[k] += chi_l * psi_Q;
		chip_U[k] += chi_l * psi_U;
		chip_V[k] += chi_l * psi_V;
	      }
	    }
	  }
	}
      }
    }
  }

  free(pf);
  return backgrflags;
}
/* ------- end ---------------------------- rlk_opacity.c ----------- */

/* ------- begin -------------------------- RLKProfile.c ------------ */

double RLKProfile(RLK_Line *rlk, int k, int mu, bool_t to_obs,
                  double lambda,
		  double *phi_Q, double *phi_U, double *phi_V,
		  double *psi_Q, double *psi_U, double *psi_V)
{
  register int nz;

  double v, phi_sm, phi_sp, phi_pi, psi_sm, psi_sp, psi_pi, adamp,
         vB, H, F, sv, phi_sigma, phi_delta, sign, sin2_gamma, phi,
         psi_sigma, psi_delta, vbroad, vtherm, GvdW, *np;
  Element *element;

  /* --- Returns the normalized profile for a Kurucz line
         and calculates the Stokes profile components if necessary -- */

  element = &atmos.elements[rlk->pt_index - 1];
  vtherm  = 2.0*KBOLTZMANN/(AMU * element->weight);
  vbroad  = sqrt(vtherm*atmos.T[k] + SQ(atmos.vturb[k]));

  v = (lambda/rlk->lambda0 - 1.0) * CLIGHT/vbroad;
  if (atmos.moving) {
    if (to_obs)
      v += vproject(k, mu) / vbroad;
    else
      v -= vproject(k, mu) / vbroad;
  }
  sv = 1.0 / (SQRTPI * vbroad);

  if (rlk->Grad) {
    switch (rlk->vdwaals) {
    case UNSOLD:
      GvdW = rlk->cross * pow(atmos.T[k], 0.3);
      break;

    case BARKLEM:
      GvdW = rlk->cross * pow(atmos.T[k], (1.0 - rlk->alpha)/2.0);
      break;

    default:
      GvdW = rlk->GvdWaals;
      break;
    }
    np = atmos.H->n[atmos.H->Nlevel-1];
    adamp = (rlk->Grad + rlk->GStark * atmos.ne[k] + 
	     GvdW * (atmos.nHtot[k] - np[k])) * 
      (rlk->lambda0  * NM_TO_M) / (4.0*PI * vbroad);
  } else {
    phi = (fabs(v) <= MAX_GAUSS_DOPPLER) ? exp(-v*v) : 0.0;
    return phi * sv;
  }

  if (rlk->polarizable) {
    sin2_gamma = 1.0 - SQ(atmos.cos_gamma[mu][k]);
    vB   = (LARMOR * rlk->lambda0) * atmos.B[k] / vbroad;
    sign = (to_obs) ? 1.0 : -1.0;

    phi_sm = phi_pi = phi_sp = 0.0;
    psi_sm = psi_pi = psi_sp = 0.0;

    for (nz = 0;  nz < rlk->zm->Ncomponent;  nz++) {
      H = Voigt(adamp, v - rlk->zm->shift[nz]*vB, &F, HUMLICEK);

      switch (rlk->zm->q[nz]) {
      case -1:
	phi_sm += rlk->zm->strength[nz] * H;
	psi_sm += rlk->zm->strength[nz] * F;
	break;
      case  0:
	phi_pi += rlk->zm->strength[nz] * H;
	psi_pi += rlk->zm->strength[nz] * F;
	break;
      case  1:
	phi_sp += rlk->zm->strength[nz] * H;
	psi_sp += rlk->zm->strength[nz] * F;
      }
    }
    phi_sigma = phi_sp + phi_sm;
    phi_delta = 0.5*phi_pi - 0.25*phi_sigma;

    phi = (phi_delta*sin2_gamma + 0.5*phi_sigma) * sv;

    *phi_Q = sign * phi_delta * sin2_gamma * atmos.cos_2chi[mu][k] * sv;
    *phi_U = phi_delta * sin2_gamma * atmos.sin_2chi[mu][k] * sv;
    *phi_V = sign * 0.5*(phi_sp - phi_sm) * atmos.cos_gamma[mu][k] * sv;

    if (input.magneto_optical) {
      psi_sigma = psi_sp + psi_sm;
      psi_delta = 0.5*psi_pi - 0.25*psi_sigma;

      *psi_Q = sign * psi_delta * sin2_gamma * atmos.cos_2chi[mu][k] * sv;
      *psi_U = psi_delta * sin2_gamma * atmos.sin_2chi[mu][k] * sv;
      *psi_V = sign * 0.5*(psi_sp - psi_sm) * atmos.cos_gamma[mu][k] * sv;
    }
  } else
   phi = Voigt(adamp, v, NULL, ARMSTRONG) * sv;

  return phi;
}
/* ------- end ---------------------------- RLKProfile.c ------------ */

/* ------- begin -------------------------- RLKZeeman.c ------------- */

void RLKZeeman(RLK_Line *rlk)
{
  const char routineName[] = "RLKZeeman";

  register int n;

  double Jl, Ju, Mu, Ml, norm[3], gLu, gLl, g_eff;

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

  Jl = rlk->level_i.J;
  Ju = rlk->level_j.J;
  
  rlk->zm = (ZeemanMultiplet *) malloc(sizeof(ZeemanMultiplet));
  initZeeman(rlk->zm);
  zm = rlk->zm;

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
  g_eff = 0.0;

  /* --- Fill the structure and normalize the strengths -- -------- */

  gLl = RLKLande(&rlk->level_i);
  gLu = RLKLande(&rlk->level_j);

  n = 0;
  for (Ml = -Jl;  Ml <= Jl;  Ml++) {
    for (Mu = -Ju;  Mu <= Ju;  Mu++) {
      if (fabs(Mu - Ml) <= 1.0) {
	zm->q[n]        = (int) (Ml - Mu);
	zm->shift[n]    = gLl*Ml - gLu*Mu;
	zm->strength[n] = ZeemanStrength(Ju, Mu, Jl, Ml);
	  
	norm[zm->q[n]+1] += zm->strength[n];
	if (zm->q[n] == 1) g_eff += zm->shift[n] * zm->strength[n];
	n++;
      }
    }
  }
  for (n = 0;  n < zm->Ncomponent;  n++)
    zm->strength[n] /= norm[zm->q[n]+1];

  g_eff /= norm[2];
  zm->g_eff = g_eff;
}
/* ------- end ---------------------------- RLKZeeman.c ------------- */

/* ------- begin -------------------------- RLKdet_level.c ---------- */

#define LS_COUNT  2
#define JK_COUNT  6
#define JJ_COUNT  4

bool_t RLKdet_level(char* label, RLK_level *level)
{
  const char routineName[] = "RLKdet_level";

  char **words, orbit[2], Lchar, L1char, lchar, Kchar, l1char, l2char;
  char delims_i[] = "({[";
  char delims_f[] = ")}]";
  char *ptr_i, *ptr_f, *quanta_str;
  bool_t counterror;
  int  count, multiplicity, length, Nread, M1;
  double J1, j1, j2;

  /* --- Get spin and orbital quantum numbers from level labels --

    For explicit LS_COUPLING, or JK_ and JJ_COUPLING provide labels with
    explicit quantum numbers as follows (note the delimiters, no spaces).
    In the Kurucz line list file (note J is always explicitly listed).

    In keyword.input set RLK_EXPLICIT = TRUE
    You CAN NOT mix explicit and non-explicit label modes!

    See:  Landi degl'Innocenti & Landolfi 2004, pp 76-77
    Note: the label HAS to be 10 characters long, and no other
          elements of the line in the .kur should be moved!!

      LS_COUPLING: [m,L]             Exmpl: '[2P]      ' --> S = 0.5, L = 1
      JK_COUPLING: (M1L1J1)lmK       Exmpl: '(6D4.5)f2K' --> S1 = 2.5, L1 = 2,
                                               J1 = 4.5, l = 3, K = 3.5
      JJ_COUPLING: {j1l1j2l2}        Exmpl: '{1.5p2.5s}' --> j1 = 1.5,
                                               l1 = 1, j2 = 2.5, l2 = 0

    Example for the FeI 1565.2874 line (JK_COUPLING in upper level):

1565.2874 -0.476 26.00   50377.905  5.0 s6D)4d f7D   56764.763  4.0 s6D9/4f[3]   8.44 -5.00 -7.70K94  0 0  0 0.000  0 0.000    0    0           1510 1542

    becomes:

1565.2874 -0.476 26.00   50377.905  5.0 [7D]         56764.763  4.0 (6D4.5)f2k   8.44 -5.00 -7.70K94  0 0  0 0.000  0 0.000    0    0           1510 1542

     --                                                -------------- */

  if (!input.RLK_explicit) {

    words = getWords(label, " ", &count);
    if (words[0]) {
      length = strlen(words[count-1]);
      Nread  = sscanf(words[count-1] + length-2, "%d%1s",
		      &multiplicity, orbit);
      if (Nread != 2 || !isupper(orbit[0])) return FALSE;

      level->L = getOrbital(orbit[0]);
      level->S = (multiplicity - 1) / 2.0;
      
      length = strlen(words[0]);
      Nread  = sscanf(words[0] + length-1, "%1s", orbit);
      free(words);

      // --- JdlCR: also extract l --- //
      level->l =  getOrbital(toupper(orbit[0]));
      
      level->cpl = LS_COUPLING;
      return TRUE;
    } else
      return FALSE;
  } else {
  
    /* --- Check for presence of any of three allowed delimiters -- - */
  
    ptr_i = strpbrk(label, delims_i);
    ptr_f = strpbrk(label, delims_f);

    if (ptr_i == NULL || ptr_f == NULL) {
      if (ptr_i != NULL && ptr_f == NULL) {
	sprintf(messageStr, " Malformed label: missing ending "
		"delimiter: %c, label: %s\n",
		delims_f[strchr(delims_i, ptr_i[0]) - delims_i], label);
	Error(ERROR_LEVEL_2, routineName, messageStr);
      }
      if (ptr_i == NULL && ptr_f != NULL) {
	sprintf(messageStr, " Malformed label: missing beginning "
		"delimiter: %c, label: %s\n",
		delims_i[strchr(delims_f, ptr_f[0]) - delims_f], label);
	Error(ERROR_LEVEL_2, routineName, messageStr);
      }
      return FALSE;
    }
    if (ptr_i[0] != delims_i[strchr(delims_f, ptr_f[0]) - delims_f]) {

      /* --- When delimiters are not matching --       -------------- */
      
      sprintf(messageStr,
	      " Malformed label: mismatched delimiters: '%c %c', %s\n",
	      ptr_i[0], ptr_f[0], label);
      Error(ERROR_LEVEL_2, routineName, messageStr);
      return FALSE;
    }     

    length = ptr_f - ptr_i + 2;
    quanta_str = (char *) malloc(length);
    strncpy(quanta_str, ptr_i, length-1);
    quanta_str[length-1] = '\0';

    counterror = FALSE;

    switch (ptr_i[0]) {
    case '[':
      if (Nread = sscanf(quanta_str, "[%1d%1c]",
			 &multiplicity, &Lchar) != LS_COUNT)
	counterror = TRUE;
      else {
	level->S = (multiplicity - 1) / 2.0;
	level->L = getOrbital(Lchar);
	
	level->cpl = LS_COUPLING;
      }
      break;

    case '(':
      if (Nread = sscanf(label, "(%1d%1c%3lf)%1c%1d%1c",
			 &M1, &L1char, &level->J1, &lchar,
			 &multiplicity, &Kchar) != JK_COUNT)
	counterror = TRUE;
      else {
	level->S1 = (M1 - 1) / 2.0;
	level->L1 = getOrbital(L1char);
	level->l  = getOrbital(toupper(lchar));
	level->K  = getJK_K(Kchar);
	
	level->cpl = JK_COUPLING;
      }
      break;

    case '{':
      if (Nread = sscanf(quanta_str, "{%3lf%1c%3lf%1c}",
			 &level->j1, &l1char,
			 &level->j2, &l2char) != JJ_COUNT)
	counterror = TRUE;
      else {
	level->l1 = getOrbital(toupper(l1char));
	level->l2 = getOrbital(toupper(l2char));
	
	level->cpl = JJ_COUPLING;
      }
      break;
    }
    if (counterror) {
      sprintf(messageStr, "Wrong quantum number count: %s: %d\n",
	      quanta_str, Nread);
      Error(ERROR_LEVEL_2, routineName, messageStr);
      return FALSE;
    }
    free(quanta_str);
    return TRUE;
  }
}
/* ------- end ---------------------------- RLKdet_level.c ---------- */

/* ------- begin -------------------------- initRLK.c --------------- */

void initRLK(RLK_Line *rlk)
{
  rlk->polarizable = FALSE;
  rlk->zm = NULL;
  rlk->level_i.l = 0; // JdlCR: init to zero
  rlk->level_j.l = 0; // JdlCR: init to zero
  rlk->cross = 0; // JdlCR: init to zero
  rlk->alpha = 0; // JdlCR: init to zero
}
/* ------- end ---------------------------- initRLK.c --------------- */

/* ------- begin -------------------------- getUnsoldcross.c -------- */

void getUnsoldcross(RLK_Line *rlk)
{
  const double FOURPIEPS0 = 4.0 * PI * EPSILON_0;

  double   Z, deltaR, vrel35_H, vrel35_He, C625;
  Element *element, *He;

  element = &atmos.elements[rlk->pt_index - 1];
  He = &atmos.elements[1];

  if (rlk->stage > element->Nstage - 1) {
    rlk->vdwaals = KURUCZ;
    return;
  }
  
  Z = rlk->stage + 1;
  deltaR = SQ(E_RYDBERG/(element->ionpot[rlk->stage] - rlk->level_j.E)) -
    SQ(E_RYDBERG/(element->ionpot[rlk->stage] - rlk->level_i.E));

  if (deltaR <= 0.0) {
    rlk->vdwaals = KURUCZ;
    return;
  }

  vrel35_H  = pow(8.0*KBOLTZMANN/(PI * AMU * element->weight) * 
		  (1.0 + element->weight/atmos.H->weight), 0.3);
  vrel35_He = pow(8.0*KBOLTZMANN/(PI * AMU * element->weight) * 
		  (1.0 + element->weight/He->weight), 0.3);

  C625 = pow(2.5 * (SQ(Q_ELECTRON)/FOURPIEPS0) *
	     (ABARH/FOURPIEPS0) *
	     2*PI * SQ(Z*RBOHR)/HPLANCK * deltaR, 0.4);
  rlk->cross = 8.08 *(vrel35_H + He->abund*vrel35_He) * C625;

  rlk->vdwaals = UNSOLD;
}
/* ------- end ---------------------------- getUnsoldcross.c -------- */

/* ------- begin -------------------------- free_BS.c --------------- */

void free_BS(Barklemstruct *bs)
{
  free(bs->neff1);
  free(bs->neff2);

  freeMatrix((void **) bs->cross);
  freeMatrix((void **) bs->alpha);
}
/* ------- end ---------------------------- free_BS.c --------------- */

/* ------- begin -------------------------- RLKLande.c -------------- */

double RLKLande(RLK_level *level)
{
  const char routineName[] = "RLKLande";

  /* --- Lande g factors for different angular momentum coupling 
         schemes.

         See: Landi degl'Innocenti & Landolfi 2004, pp 76-77 -- ----- */

  switch (level->cpl) {
  case LS_COUPLING:
    level->gL = 1.0 + zm_gamma(level->J, level->S, level->L);

    break;
    
  case JK_COUPLING:
    level->gL = 1.0 + zm_gamma(level->J, 0.5, level->K) +
      zm_gamma(level->J, level->K, 0.5) *
      zm_gamma(level->K, level->J1, level->l) *
      zm_gamma(level->J1, level->S1, level->L1);

    break;
    
  case JJ_COUPLING:
    level->gL = 1.0 + zm_gamma(level->J, level->j1, level->j2) *
      zm_gamma(level->j1, 0.5, level->l1) +
      zm_gamma(level->J, level->j2, level->j1) *
      zm_gamma(level->j2, 0.5, level->l2);

    break;
    
  default:
    sprintf(messageStr, "Invalid coupling: %d", level->cpl);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }
  return level->gL;
}
/* ------- end ---------------------------- RLKLande.c -------------- */

/* ------- begin -------------------------- getJK_K.c -- ------------ */

double getJK_K(char Kchar)
{
  const char routineName[] = "getJK_K";

  double K;

  switch (Kchar) {
  case 'p': K = 0.5;  break;
  case 'f': K = 1.5;  break;
  case 'h': K = 2.5;  break;
  case 'k': K = 3.5;  break;
  case 'm': K = 4.5;  break;
  case 'o': K = 5.5;  break;
  case 'r': K = 6.5;  break;
  case 't': K = 7.5;  break;
  case 'u': K = 8.5;  break;
  case 'v': K = 9.5;  break;
  case 'w': K = 10.5; break;
  default: 
    sprintf(messageStr, "Invalid Kchar: %c", Kchar);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }
  return K;
}
/* ------- end ---------------------------- getJK_K.c --------------- */
