/* ------- file: -------------------------- readatmos.c -------------

       Version:       rh2.0, 1-D spherically symmetric
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Mon Sep 22 16:39:01 2025 --

       --------------------------                      ----------RH-- */

/* --- Reads atmospheric model in spherical MULTI format. --  ------- */


#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "geometry.h"
#include "spectrum.h"
#include "background.h"
#include "constant.h"
#include "error.h"
#include "inputs.h"
#include "statistics.h"
#include "xdr.h"


#define MULTI_COMMENT_CHAR  "*"
#define N_HYDROGEN_MULTI     6


/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

extern Spectrum spectrum;
extern InputData input;
extern char messageStr[];


/* ------- begin -------------------------- readAtmos.c ------------- */

void readAtmos(Atmosphere *atmos, Geometry *geometry)
{
  const char routineName[] = "readAtmos";
  register int k, n;

  char    scaleStr[21], inputLine[MAX_LINE_SIZE], *filename;
  bool_t  exit_on_EOF, enhanced_atmos_ID = FALSE;
  int     Nread, Nradius, Nrequired, checkPoint;
  double  *dscale, *r, *cmass, *tau_ref, turbpress, nbaryon;
  struct  stat statBuffer;

  getCPU(2, TIME_START, NULL);

  /* --- Get abundances of background elements --        ------------ */
 
  readAbundance(atmos);

  /* --- Open input file for model atmosphere in SMULTI format -- --- */

  if ((atmos->fp_atmos = fopen(input.atmos_input, "r")) == NULL) {
    sprintf(messageStr, "Unable to open inputfile %s",
	    input.atmos_input);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  } else {
    sprintf(messageStr, " -- reading input file: %s\n", input.atmos_input);
    Error(MESSAGE, NULL, messageStr);
  }

  atmos->NHydr = N_HYDROGEN_MULTI;

  /* --- Boundary condition at TOP of atmosphere --      ------------ */

  if (strcmp(input.Itop, "none"))
    geometry->rboundary[TOP] = IRRADIATED;
  else 
    geometry->rboundary[TOP] = ZERO;

  /* --- Boundary condition at BOTTOM of atmosphere --   ------------ */

  geometry->rboundary[CORE] = THERMALIZED;

  /* --- Read atmos ID, scale type, gravity, and number of depth
         points --                                       ------------ */
 
  getLine(atmos->fp_atmos, MULTI_COMMENT_CHAR, inputLine, exit_on_EOF=TRUE);
  if (enhanced_atmos_ID) {

    /* --- Construct atmosID from filename and last modification date */

    stat(input.atmos_input, &statBuffer);
    if ((filename = strrchr(input.atmos_input, '/')))
      filename++;
    else
      filename = input.atmos_input;
    sprintf(atmos->ID, "%s (%.24s)", filename,
	    asctime(localtime(&statBuffer.st_mtime)));
    Nread = 1;
  } else
    Nread = sscanf(inputLine, "%s", atmos->ID);

  getLine(atmos->fp_atmos, MULTI_COMMENT_CHAR, inputLine, exit_on_EOF=TRUE);
  Nread += sscanf(inputLine, "%20s", scaleStr);
  getLine(atmos->fp_atmos, MULTI_COMMENT_CHAR, inputLine, exit_on_EOF=TRUE);
  Nread += sscanf(inputLine, "%lf %lf", &atmos->gravity, &geometry->Radius);
  getLine(atmos->fp_atmos, MULTI_COMMENT_CHAR, inputLine, exit_on_EOF=TRUE);
  Nread += sscanf(inputLine, "%d %d %d", &Nradius,
		  &geometry->Ncore, &geometry->Ninter);
  checkNread(Nread, Nrequired=7, routineName, checkPoint=1);

  /* --- Keep duplicates of some of the geometrical quantities in
         Atmos structure --                            -------------- */

  atmos->Ndim = 1;
  atmos->N = (int *) malloc(atmos->Ndim * sizeof(int));
  atmos->N[0] = geometry->Nradius = Nradius;
  atmos->Nspace = Nradius;


  atmos->gravity = POW10(atmos->gravity) * CM_TO_M;
  geometry->Radius *= KM_TO_M;

  /* --- Allocate space for arrays that define structure -- --------- */

  tau_ref = geometry->tau_ref = (double *) malloc(Nradius * sizeof(double));
  cmass  = geometry->cmass  = (double *) malloc(Nradius * sizeof(double));
  r = geometry->r = (double *) malloc(Nradius * sizeof(double));
  atmos->T      = (double *) malloc(Nradius * sizeof(double));
  atmos->ne     = (double *) malloc(Nradius * sizeof(double));
  atmos->vturb  = (double *) malloc(Nradius * sizeof(double));
  geometry->vel = (double *) malloc(Nradius * sizeof(double));


  dscale = (double *) malloc(Nradius * sizeof(double));
  for (k = 0;  k < Nradius;  k++) {
    getLine(atmos->fp_atmos, MULTI_COMMENT_CHAR, inputLine, exit_on_EOF=TRUE);
    Nread = sscanf(inputLine, "%lf %lf %lf %lf %lf",
		    &dscale[k], &atmos->T[k], &atmos->ne[k],
		    &geometry->vel[k], &atmos->vturb[k]);
    checkNread(Nread, Nrequired=5, routineName, checkPoint=2);
  }

  switch(toupper(scaleStr[0])) {
  case 'M':
    geometry->scale = COLUMN_MASS;
    for (k = 0;  k < Nradius;  k++)
      geometry->cmass[k] = POW10(dscale[k]) * (G_TO_KG / SQ(CM_TO_M));
    break;
  case 'T':
    geometry->scale = TAU500;
    for (k = 0;  k < Nradius;  k++) geometry->tau_ref[k] = POW10(dscale[k]);
    break;
  case 'H':
    geometry->scale  = GEOMETRIC;
    for (k = 0;  k < Nradius;  k++) geometry->r[k] = dscale[k] * KM_TO_M;
    break;
  default:
    sprintf(messageStr, "Unknown depth scale string in file %s: %s",
	    input.atmos_input, scaleStr);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }
  free(dscale);

  for (k = 0;  k < Nradius;  k++) {
    geometry->vel[k] *= KM_TO_M;
    atmos->vturb[k]  *= KM_TO_M;
    atmos->ne[k]     /= CUBE(CM_TO_M);
  }
  atmos->moving = FALSE;
  for (k = 0;  k < Nradius;  k++) {
    if (fabs(geometry->vel[k]) > atmos->vmacro_tresh) {
      atmos->moving = TRUE;
      break;
    }
  }
  if (atmos->moving)
    Error(ERROR_LEVEL_2, routineName,
	  "Cannot yet deal with macroscopic flows in spherical geometry");
    
  /* --- Read Hydrogen populations --                     ----------- */

  atmos->nH = matrix_double(atmos->NHydr, Nradius);
  for (k = 0;  k < Nradius;  k++) {
    if (getLine(atmos->fp_atmos, MULTI_COMMENT_CHAR, inputLine,
		exit_on_EOF=FALSE) == EOF) break;
    Nread = sscanf(inputLine, "%lf %lf %lf %lf %lf %lf",
		   &atmos->nH[0][k], &atmos->nH[1][k], &atmos->nH[2][k],
		   &atmos->nH[3][k], &atmos->nH[4][k], &atmos->nH[5][k]);
    checkNread(Nread, Nrequired=6, routineName, checkPoint=3);
  }

  if (k > 0  &&  k < Nradius) {
    sprintf(messageStr,
	    "Reached end of input file %s before all data was read",
	    input.atmos_input);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  } else if (k == 0) {

    /* --- No hydrogen populations supplied: use LTE populations
           like MULTI does --                          -------------- */

    if (geometry->scale != COLUMN_MASS) {
      sprintf(messageStr,
	      "Height scale should be COLUMNMASS when nH not supplied: "
	      "File %s", input.atmos_input);
      Error(ERROR_LEVEL_2, routineName, messageStr);
    }
    atmos->nHtot = (double *) calloc(Nradius, sizeof(double));

    atmos->H_LTE = TRUE;    
    for (k = 0;  k < Nradius;  k++) {
      turbpress = 0.5 * atmos->avgMolWght*AMU * SQ(atmos->vturb[k]);
      nbaryon =
	(atmos->gravity * geometry->cmass[k] + turbpress * atmos->ne[k]) /
	(KBOLTZMANN * atmos->T[k] + turbpress) - atmos->ne[k];
      atmos->nHtot[k] = nbaryon / atmos->totalAbund;
    }
  } else if (k == Nradius) {
    atmos->nHtot = (double *) calloc(Nradius, sizeof(double));
    for (n = 0;  n < atmos->NHydr;  n++) {
      for (k = 0;  k < Nradius;  k++) {
	atmos->nH[n][k] /= CUBE(CM_TO_M);
	atmos->nHtot[k] += atmos->nH[n][k];
      }
    }
  }
  /* --- Magnetic field is read here. --               -------------- */

  atmos->Stokes = FALSE;

  getCPU(2, TIME_POLL, "Read Atmosphere");
}
/* ------- end ---------------------------- readAtmos.c ------------- */

/* ------- begin -------------------------- convertScales.c --------- */

void convertScales(Atmosphere *atmos, Geometry *geometry)
{
  register int k;

  bool_t hunt;
  int    ref_index, Nradius = geometry->Nradius;
  double *rho, *r, *cmass, *tau_ref, r_zero, unity;
  ActiveSet *as;

  /* --- Convert between different depth scales --       ------------ */

  r       = geometry->r;
  cmass   = geometry->cmass;
  tau_ref = geometry->tau_ref;

  rho = (double *) malloc(Nradius * sizeof(double));
  for (k = 0;  k < Nradius;  k++)
    rho[k] = (AMU * atmos->wght_per_H) * atmos->nHtot[k];

  /* --- Get opacity of reference wavelength --          ------------ */

  Locate(spectrum.Nspect, spectrum.lambda, atmos->lambda_ref, &ref_index);

  as = &spectrum.as[ref_index];
  alloc_as(ref_index, FALSE);
  readBackground(ref_index, 0, FALSE);

  /* --- Convert to missing depth scales --              ------------ */

  switch (geometry->scale) {
  case COLUMN_MASS:
    r[0] = 0.0;
    tau_ref[0] = as->chi_c[0] / rho[0] * cmass[0];
    for (k = 1;  k < Nradius;  k++) {
      r[k] = r[k-1] - 2.0*(cmass[k] - cmass[k-1]) / 
	(rho[k-1] + rho[k]);
      tau_ref[k] = tau_ref[k-1] + 0.5*(as->chi_c[k-1] + as->chi_c[k]) *
	(r[k-1] - r[k]);
    }
    break;
  case TAU500:
    r[0] = 0.0;
    cmass[0]  = (tau_ref[0] / as->chi_c[0]) * rho[0];
    for (k = 1;  k < Nradius;  k++) {
      r[k] = r[k-1] - 2.0 * (tau_ref[k] - tau_ref[k-1]) / 
	(as->chi_c[k-1] + as->chi_c[k]);
      cmass[k]  = cmass[k-1]  + 0.5*(rho[k-1] + rho[k]) *
	(r[k-1] - r[k]);
    }
    break;
  case GEOMETRIC:
    cmass[0] = tau_ref[0] = 0.0;
    for (k = 1;  k < Nradius;  k++) {
      cmass[k]  = cmass[k-1]  + 0.5*(rho[k-1] + rho[k]) *
	(r[k-1] - r[k]);
      tau_ref[k] = tau_ref[k-1] + 0.5*(as->chi_c[k-1] + as->chi_c[k]) *
	(r[k-1] - r[k]);
    }
    break;
  }
  free_as(ref_index, FALSE);

  if ((geometry->scale == COLUMN_MASS) || (geometry->scale == TAU500)) {
    unity = 1.0;
    Linear(Nradius, tau_ref, r, 1, &unity, &r_zero, hunt=FALSE);
    for (k = 0;  k < Nradius;  k++)  r[k] = r[k] - r_zero;
  }
  free(rho);
}
/* ------- end ---------------------------- convertScales.c --------- */

/* ------- begin -------------------------- getBoundary.c ----------- */

void getBoundary(Geometry *geometry)
{
  const char routineName[] = "getBoundary";
  register int la;

  bool_t result = TRUE;
  FILE  *fp_Itop;
  XDR    xdrs;

  switch (geometry->rboundary[TOP]) {
  case ZERO: break;
  case IRRADIATED:
    sprintf(messageStr, "\n -- reading irradiance input file: %s\n\n",
	    input.Itop);
    Error(MESSAGE, NULL, messageStr);

    geometry->Itop = (double *) malloc(spectrum.Nspect * sizeof(double));

    /* --- Open input file for irradiation at TOP --     -------------- */
    
    if ((fp_Itop = fopen(input.Itop, "r")) == NULL) {
      sprintf(messageStr, "Unable to open inputfile %s", input.Itop);
      Error(ERROR_LEVEL_2, routineName, messageStr);
    }
    xdrstdio_create(&xdrs, fp_Itop, XDR_DECODE);

    result &= xdr_vector(&xdrs, (char *) geometry->Itop,
			 spectrum.Nspect, sizeof(double),
			 (xdrproc_t) xdr_double);
    if (!result) {
      sprintf(messageStr,
	      "Unable to read irradiation data at TOP of atmosphere");
      Error(ERROR_LEVEL_2, routineName, messageStr);
    }
    xdr_destroy(&xdrs);
    fclose(fp_Itop);
    break;
  default:
    Error(ERROR_LEVEL_2, routineName,
	  "Invalid boundary condition at the TOP of atmosphere");
  }

  switch (geometry->rboundary[CORE]) {
  case THERMALIZED:
    break;
  default:
    Error(ERROR_LEVEL_2, routineName,
	  "Invalid boundary condition at the BOTTOM of atmosphere");
  }

}
/* ------- end ---------------------------- getBoundary.c ----------- */
