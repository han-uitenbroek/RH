/* ------- file: -------------------------- writegeom_xdr.c ---------

       Version:       rh2.0, 1-D spherically symmetric
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Fri Jun  2 11:34:00 2000 --

       --------------------------                      ----------RH-- */

/* --- Writes geometry data for 1-D spherical-symmetric version
       to output.

       XDR (external data representation) version. --  -------------- */

 
#include <stdlib.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "geometry.h"
#include "error.h"
#include "inputs.h"
#include "xdr.h"

/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

extern enum Topology topology;

extern InputData  input;
extern char   messageStr[];


/* ------- begin -------------------------- writeGeometry.c --------- */

void writeGeometry(Geometry *geometry)
{
  const char routineName[] = "writeGeometry";
  register int mu;

  bool_t  result;
  int     Nrays = geometry->Nrays, Nradius = geometry->Nradius;
  FILE   *fp_out;
  double *xmu0, *wmu0;
  XDR     xdrs;

  if (!strcmp(input.geometry_output, "none")) return;

  if ((fp_out = fopen(input.geometry_output, "w")) == NULL) {
    sprintf(messageStr, "Unable to open output file %s",
	    input.geometry_output);
    Error(ERROR_LEVEL_1, routineName, messageStr);
    return;
  }
  xdrstdio_create(&xdrs, fp_out, XDR_ENCODE);

  result = xdr_enum(&xdrs, (enum_t *) &topology);
  result = xdr_int(&xdrs, &geometry->Nrays);
  result = xdr_int(&xdrs, &geometry->Nradius);
  result = xdr_int(&xdrs, &geometry->Ncore);

  result = xdr_double(&xdrs, &geometry->Radius);
  xmu0 = (double *) malloc(Nrays * sizeof(double));
  wmu0 = (double *) malloc(Nrays * sizeof(double));
  for (mu = 0;  mu < Nrays;  mu++) {
    xmu0[mu] = geometry->rays[mu].xmu[0];
    wmu0[mu] = geometry->rays[mu].wmu[0];
  }
  result = xdr_vector(&xdrs, (char *) xmu0, Nrays,
		      sizeof(double), (xdrproc_t) xdr_double);
  result = xdr_vector(&xdrs, (char *) wmu0, Nrays,
		      sizeof(double), (xdrproc_t) xdr_double);
  result = xdr_vector(&xdrs, (char *) geometry->r, Nradius,
		      sizeof(double), (xdrproc_t) xdr_double);
  result = xdr_vector(&xdrs, (char *) geometry->cmass, Nradius,
		      sizeof(double), (xdrproc_t) xdr_double);
  result = xdr_vector(&xdrs, (char *) geometry->tau_ref, Nradius,
		      sizeof(double), (xdrproc_t) xdr_double);
  result = xdr_vector(&xdrs, (char *) geometry->vel, Nradius,
		       sizeof(double),  (xdrproc_t) xdr_double);

  xdr_destroy(&xdrs);
  fclose(fp_out);
  free(xmu0);  free(wmu0);
}
/* ------- end ---------------------------- writeGeometry.c --------- */
