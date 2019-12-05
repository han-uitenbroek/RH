/* ------- file: -------------------------- writegeom_xdr.c ---------

       Version:       rh2.0, 3-D Cartesian
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Dec 13 1996

       --------------------------                      ----------RH-- */

/* --- Writes geometry data for 3-D Cartesian version to output file.
       XDR (external data representation) version. --  -------------- */

 
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

extern Atmosphere atmos;
extern InputData input;
extern char messageStr[];


/* ------- begin -------------------------- writeGeometry.c --------- */

void writeGeometry(Geometry *geometry)
{
  const char routineName[] = "writeGeometry";

  int   result, Nrays = geometry->Nrays,
        Nspace = (geometry->Nx * geometry->Ny * geometry->Nz);
  FILE *fp_out;
  XDR   xdrs;

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
  result = xdr_int(&xdrs, &geometry->Nx);
  result = xdr_int(&xdrs, &geometry->Ny);
  result = xdr_int(&xdrs, &geometry->Nz);
  result = xdr_enum(&xdrs, (enum_t *) &atmos.angleSet.set);

  result = xdr_vector(&xdrs, (char *) geometry->mux, Nrays,
		      sizeof(double),  (xdrproc_t) xdr_double);
  result = xdr_vector(&xdrs, (char *) geometry->muy, Nrays,
		      sizeof(double), (xdrproc_t) xdr_double);
  result = xdr_vector(&xdrs, (char *) geometry->wmu, Nrays,
		      sizeof(double), (xdrproc_t) xdr_double);

  result = xdr_double(&xdrs, &geometry->dx);
  result = xdr_double(&xdrs, &geometry->dy);
  result = xdr_vector(&xdrs, (char *) geometry->z, geometry->Nz,
		      sizeof(double), (xdrproc_t) xdr_double);

  result = xdr_vector(&xdrs, (char *) geometry->vx, Nspace,
		      sizeof(double),  (xdrproc_t) xdr_double);
  result = xdr_vector(&xdrs, (char *) geometry->vy, Nspace,
		      sizeof(double),  (xdrproc_t) xdr_double);
  result = xdr_vector(&xdrs, (char *) geometry->vz, Nspace,
		      sizeof(double),  (xdrproc_t) xdr_double);

  xdr_destroy(&xdrs);
  fclose(fp_out);
}
/* ------- end ---------------------------- writeGeometry.c --------- */
