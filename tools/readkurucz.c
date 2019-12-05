/* ------- file: -------------------------- readkurucz.c ------------

       Version:       rh1.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Fri Apr 29 03:56:49 2011 --

       --------------------------                      ----------RH-- */

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atomweights.h"
#include "atmos.h"
#include "rhf1d/geometry.h"
#include "error.h"
#include "inputs.h"
#include "constant.h"


#define COMMENT_CHAR  "*"


/* --- Function prototypes --                          -------------- */

void writeMULTIatmos(Geometry *geometry, Atmosphere *atmos,
		     char *modelName);


/* --- Global variables --                             -------------- */

CommandLine commandline;
char messageStr[MAX_MESSAGE_LENGTH];


/* ------- begin -------------------------- readKurucz.c ------------ */

int main(int argc, char *argv[])
{
  register int n, k;

  char   line[MAX_LINE_SIZE], *sread;
  int    Nread, no;
  double Teff, abundscale, logabundH, pressure, Rosseland_opac, rad_acc;

  Atmosphere atmos;
  Geometry geometry;
  Element *element;
  FILE  *fp_Kurucz;

  commandline.logfile = stderr;

  if (argc < 2) {
    fprintf(stderr, "Usage : %s  model.Kurucz [MULTI_name]\n\n",
	    (char *) argv[0]);
    exit(0);
  }

  atmos.Nelem    = sizeof(atomweight) / sizeof(struct AtomWeight);
  atmos.elements = (Element *) malloc(atmos.Nelem * sizeof(Element));

  for (n = 0;  n < atmos.Nelem;  n++) {
    element = &atmos.elements[n];
    strcpy(element->ID, atomweight[n].ID);
    element->weight = atomweight[n].weight;
  }

  if ((fp_Kurucz = fopen(argv[1], "r")) == NULL) {
    sprintf(messageStr, "Unable to open input file %s",
	    (char *) argv[1]);
    Error(ERROR_LEVEL_2, argv[0], messageStr);
  }

  sread = fgets(line, MAX_LINE_SIZE, fp_Kurucz);
  Nread = sscanf(line, " TEFF %lf GRAVITY %lf", &Teff, &atmos.gravity);
  atmos.gravity = POW10(atmos.gravity);

  sread = fgets(line, MAX_LINE_SIZE, fp_Kurucz);
  strncpy(atmos.ID, line+5, ATMOS_ID_WIDTH);
  atmos.ID[ATMOS_ID_WIDTH-1] = '\0';

  for (n = 0;  n < 3;  n++)
    sread = fgets(line, MAX_LINE_SIZE, fp_Kurucz);
  Nread = sscanf(line,
		 "ABUNDANCE SCALE %lf ABUNDANCE CHANGE %d %lf %d %lf",
		 &abundscale, &no, &atmos.elements[0].abund,
		 &no, &atmos.elements[1].abund);

  if (abundscale != 1.0) {
    sprintf(messageStr, "Use METALICITY = %f in keyword.input\n",
	    log10(abundscale));
    Error(WARNING, argv[0], messageStr);
  }

  for (n = 2;  n < 98;  n += 6) {
    sread = fgets(line, MAX_LINE_SIZE, fp_Kurucz);
    Nread = sscanf(line,
		   "ABUNDANCE CHANGE %d %lf %d %lf %d %lf %d %lf %d %lf %d %lf",
		   &no, &atmos.elements[n].abund,
		   &no, &atmos.elements[n+1].abund,
		   &no, &atmos.elements[n+2].abund,
		   &no, &atmos.elements[n+3].abund,
		   &no, &atmos.elements[n+4].abund,
		   &no, &atmos.elements[n+5].abund);
  }
  sread = fgets(line, MAX_LINE_SIZE, fp_Kurucz);
  Nread = sscanf(line, "ABUNDANCE CHANGE %d %lf",
		 &no, &atmos.elements[98].abund);

  logabundH = log10(atmos.elements[0].abund);
  atmos.elements[0].abund = 12.0;
  atmos.elements[1].abund = log10(atmos.elements[1].abund) - logabundH;
  for (n = 2;  n < atmos.Nelem;  n++)
    atmos.elements[n].abund -= logabundH;

  sread = fgets(line, MAX_LINE_SIZE, fp_Kurucz);
  Nread = sscanf(line, "READ DECK6 %d", &geometry.Ndep);
  atmos.Nspace = geometry.Ndep;

  geometry.scale = COLUMN_MASS;
  geometry.cmass = (double *) malloc(geometry.Ndep * sizeof(double));
  atmos.T        = (double *) malloc(geometry.Ndep * sizeof(double));
  atmos.ne       = (double *) malloc(geometry.Ndep * sizeof(double));
  geometry.vel   = (double *) calloc(geometry.Ndep, sizeof(double));
  atmos.vturb    = (double *) malloc(geometry.Ndep * sizeof(double));

  for (k = 0;  k < geometry.Ndep;  k++) {
    sread = fgets(line, MAX_LINE_SIZE, fp_Kurucz);
    Nread = sscanf(line, "%lf %lf %lf %lf %lf %lf %lf",
		   &geometry.cmass[k], &atmos.T[k], &pressure,
		   &atmos.ne[k], &Rosseland_opac, &rad_acc, &atmos.vturb[k]);

    atmos.vturb[k] *= CM_TO_M;
  }
  atmos.H_LTE = TRUE;

  atmos.B = NULL;

  fclose(fp_Kurucz);
  if (argc == 1) 
    writeMULTIatmos(&geometry, &atmos, NULL);
  else
    writeMULTIatmos(&geometry, &atmos, argv[2]);
}
/* ------- end ---------------------------- readKurucz.c ------------ */
