/* ------- file: -------------------------- writemulti.c ------------

       Version:       rh1.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Thu Oct  1 16:48:29 2009 --

       --------------------------                      ----------RH-- */

/* --- Write atmospheric model in MULTI format given the structures
       atmos and geometry --                           -------------- */

#include <math.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "rhf1d/geometry.h"
#include "error.h"
#include "constant.h"
#include "xdr.h"


#define COMMENT_CHAR  "*"

/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */


/* ------- begin -------------------------- writeMULTIatmos.c ------- */

void writeMULTIatmos(Geometry *geometry, Atmosphere *atmos,
		     char *modelName)
{
  register int k;

  char    fileName[MAX_LINE_SIZE];
  bool_t  result = TRUE;
  double *heightscale, height, vturb, *np;
  FILE   *fp_MULTI;
  XDR     xdrs;

  if (modelName != NULL) {
    sprintf(fileName, "%s%s", modelName, ".atmos");
    fp_MULTI = fopen(fileName, "w");
    fprintf(fp_MULTI, "%s %s%s\n  %s\n", COMMENT_CHAR, atmos->ID,
	    COMMENT_CHAR, modelName);
  } else
    fp_MULTI = stdout;

  switch(geometry->scale) {
  case COLUMN_MASS:
    fprintf(fp_MULTI, "  Mass scale\n");
    break; 
  case TAU500:
    fprintf(fp_MULTI, "  Tau 500\n");
    break;
  case GEOMETRIC:
    fprintf(fp_MULTI, "  Height\n");
  }

  fprintf(fp_MULTI, "%s\n%s log g [cm s^-2]\n %7.4f\n",
	  COMMENT_CHAR, COMMENT_CHAR, log10(atmos->gravity));
  fprintf(fp_MULTI, "%s\n%s Ndep\n %3d\n",
	  COMMENT_CHAR, COMMENT_CHAR, geometry->Ndep);

  switch(geometry->scale) {
  case COLUMN_MASS:
    fprintf(fp_MULTI, "%slg column Mass     Temperature        Ne         V"
	  "              Vturb\n", COMMENT_CHAR);
    heightscale = geometry->cmass;
    break;
  case TAU500:
    fprintf(fp_MULTI, "%slog tau 500        Temperature        Ne         V"
	  "              Vturb\n", COMMENT_CHAR);
    heightscale = geometry->tau_ref;
    break;
  case GEOMETRIC:
    fprintf(fp_MULTI, "%sheight [km]        Temperature        Ne         V"
	  "              Vturb\n", COMMENT_CHAR);
    heightscale = geometry->height;
  }

  for (k = 0;  k < geometry->Ndep;  k++) {
    if (geometry->scale == COLUMN_MASS || geometry->scale == TAU500)
      height = log10(heightscale[k]);
    vturb = atmos->vturb[k] / KM_TO_M;
    fprintf(fp_MULTI, "%17.8E %14.6E %14.6E %14.6E %14.6E\n",
	    height, atmos->T[k], atmos->ne[k], geometry->vel[k], vturb);
  }
  fprintf(fp_MULTI, "%s\n%s Hydrogen populations\n", COMMENT_CHAR, COMMENT_CHAR);
  if (atmos->H_LTE) {
    fprintf(fp_MULTI, "%s LTE values\n", COMMENT_CHAR);
  } else {
    np = atmos->nH[atmos->NHydr-1];

    fprintf(fp_MULTI,	    "%s     nh(1)       nh(2)       nh(3)"
            "       nh(4)       nh(5)       np\n", COMMENT_CHAR);
    for (k = 0;  k < geometry->Ndep;  k++) {
      fprintf(fp_MULTI, " %11.4E %11.4E %11.4E %11.4E %11.4E %11.4E\n",
	      atmos->nH[0][k], atmos->nH[1][k], atmos->nH[2][k],
	      atmos->nH[3][k], atmos->nH[4][k], np[k]);
    }
  }
  fclose(fp_MULTI);

  if (atmos->B != NULL) {
    sprintf(fileName, "%s%s", modelName, ".B");
    fp_MULTI = fopen(fileName, "w");
    xdrstdio_create(&xdrs, fp_MULTI, XDR_ENCODE);

    result &= xdr_vector(&xdrs, (char *) atmos->B, atmos->Nspace,
			 sizeof(double), (xdrproc_t) xdr_double);
    result &= xdr_vector(&xdrs, (char *) atmos->gamma_B, atmos->Nspace,
			 sizeof(double), (xdrproc_t) xdr_double);
    result &= xdr_vector(&xdrs, (char *) atmos->chi_B, atmos->Nspace,
			 sizeof(double), (xdrproc_t) xdr_double);

    xdr_destroy(&xdrs);
    fclose(fp_MULTI);
  }
}
/* ------- end ---------------------------- writeMULTIatmos.c ------- */
