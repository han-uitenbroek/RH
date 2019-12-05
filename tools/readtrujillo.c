/* ------- file: -------------------------- readtrujillo.c ----------

       Version:       rh1.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Tue Apr  4 17:30:52 2000 --

       --------------------------                      ----------RH-- */

#include <malloc.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "rhf1d/geometry.h"
#include "error.h"
#include "inputs.h"
#include "constant.h"


#define COMMENT_CHAR      "#"
#define DEGREE_TO_RADIAN  0.017453293
#define GAUSS_TO_TESLA    1.0E-4


/* --- Function prototypes --                          -------------- */

void writeMULTIatmos(Geometry *geometry, Atmosphere *atmos,
		     char *modelName);


/* --- Global variables --                             -------------- */

CommandLine commandline;
char messageStr[MAX_MESSAGE_LENGTH];


/* ------- begin -------------------------- readTrujillo.c ---------- */

void main(int argc, void *argv[])
{
  register int n, k;

  char   line[MAX_LINE_SIZE];
  bool_t exit_on_EOF;
  int    Nread;
  double Pe;

  Atmosphere atmos;
  Geometry geometry;
  FILE  *fp_Trujillo, *fp_MULTI;

  commandline.logfile = stdout;
  if (argc < 3) {
    fprintf(stderr, "Usage : %s  Trujillo.model  MULTI_name\n\n", argv[0]);
    exit(0);
  }

  if ((fp_Trujillo = fopen(argv[1], "r")) == NULL) {
    sprintf(messageStr, "Unable to open input file %s", argv[1]);
    Error(ERROR_LEVEL_2, argv[0], argv[1]);
  }

  strncpy(atmos.ID, "Atmospheric model from Javier Trujillo Bueno",
	  ATMOS_ID_WIDTH);
  atmos.ID[ATMOS_ID_WIDTH-1] = '\0';

  getLine(fp_Trujillo, COMMENT_CHAR, line, exit_on_EOF=TRUE);
  Nread = sscanf(line, "%d", &geometry.Ndep);
  atmos.Nspace = geometry.Ndep;

  atmos.gravity    = POW10(4.44);
  geometry.scale   = TAU500;
  geometry.tau_ref = (double *) malloc(geometry.Ndep * sizeof(double));

  atmos.T      = (double *) malloc(geometry.Ndep * sizeof(double));
  atmos.ne     = (double *) malloc(geometry.Ndep * sizeof(double));
  geometry.vel = (double *) malloc(geometry.Ndep * sizeof(double));
  atmos.vturb  = (double *) malloc(geometry.Ndep * sizeof(double));

  atmos.B       = (double *) malloc(geometry.Ndep * sizeof(double));
  atmos.gamma_B = (double *) malloc(geometry.Ndep * sizeof(double));
  atmos.chi_B   = (double *) malloc(geometry.Ndep * sizeof(double));

  for (k = geometry.Ndep-1;  k >= 0;  k--) {
    getLine(fp_Trujillo, COMMENT_CHAR, line, exit_on_EOF=TRUE);
    Nread = sscanf(line, "%lf %lf %lf %lf %lf %lf %lf %lf",
		   &geometry.tau_ref[k], &atmos.T[k], &Pe, &atmos.vturb[k],
		   &atmos.B[k], &geometry.vel[k], &atmos.gamma_B[k],
		   &atmos.chi_B[k]);

    geometry.tau_ref[k] = POW10(geometry.tau_ref[k]);
    atmos.ne[k] = (Pe * ERG_TO_JOULE) / (KBOLTZMANN * atmos.T[k]);

    geometry.vel[k] *= CM_TO_M;
    atmos.vturb[k]  *= CM_TO_M;

    atmos.B[k]       *= GAUSS_TO_TESLA;
    atmos.gamma_B[k] *= DEGREE_TO_RADIAN;
    atmos.chi_B[k]   *= DEGREE_TO_RADIAN;
  }
  atmos.H_LTE = TRUE;

  fclose(fp_Trujillo);
  writeMULTIatmos(&geometry, &atmos, argv[2]);
}
/* ------- end ---------------------------- readTrujillo.c ---------- */
