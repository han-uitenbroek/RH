/* ------- file: -------------------------- writeatmos_xdr.c --------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Tue Feb  8 10:57:39 2000 --

       --------------------------                      ----------RH-- */

/* --- Writes atmospheric data to output file.
       XDR (external data representation) version. --  -------------- */

 
#include <stdlib.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "error.h"
#include "inputs.h"
#include "xdr.h"

/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

extern InputData input;
extern char messageStr[];


/* ------- begin -------------------------- xdr_counted_string.c ---- */

bool_t xdr_counted_string(XDR *xdrs, char **p)
{
  bool_t output = (xdrs->x_op == XDR_ENCODE) ? TRUE : FALSE; 
  short  length;

  /* --- Reads or writes string of length length from or to xdr stream
         pointed to by xdrs.

    See: IDL User's Guide 17-34
         --                                            -------------- */

  if (output) length = strlen(*p);

  if (!xdr_short(xdrs, &length))  return FALSE;
  if (!output) {
    *p = (char *) malloc((length + 1) * sizeof(char));
    (*p)[length] = '\0';
  }
  return (length ? xdr_string(xdrs, p, length) : TRUE);
}
/* ------- end ---------------------------- xdr_counted_string.c ---- */

/* ------- begin -------------------------- writeAtmos.c ------------ */

void writeAtmos(Atmosphere *atmos)
{
  const char routineName[] = "writeAtmos";
  register int n;

  bool_t result = TRUE;
  int    Nspace = atmos->Nspace;
  char  *elemID, *atmosID;
  FILE  *fp_out;
  XDR    xdrs;

  if (!strcmp(input.atmos_output, "none")) return;

  if ((fp_out = fopen(input.atmos_output, "w")) == NULL) {
    sprintf(messageStr, "Unable to open output file %s", input.atmos_output);
    Error(ERROR_LEVEL_1, routineName, messageStr);
    return;
  }
  xdrstdio_create(&xdrs, fp_out, XDR_ENCODE);

  result &= xdr_int(&xdrs, &atmos->H->Nlevel);
  result &= xdr_int(&xdrs, &atmos->Nelem);
  result &= xdr_bool(&xdrs, &atmos->moving);

  result &= xdr_vector(&xdrs, (char *) atmos->T, Nspace,
		       sizeof(double), (xdrproc_t) xdr_double);
  result &= xdr_vector(&xdrs, (char *) atmos->ne, Nspace,
		       sizeof(double), (xdrproc_t) xdr_double);
  result &= xdr_vector(&xdrs, (char *) atmos->vturb, Nspace,
		       sizeof(double), (xdrproc_t) xdr_double);
  result &= xdr_vector(&xdrs, (char *) atmos->H->n[0],
		       Nspace * atmos->H->Nlevel,
		       sizeof(double), (xdrproc_t) xdr_double);

  atmosID = atmos->ID;
  result &= xdr_counted_string(&xdrs, &atmosID);

  for (n = 0;  n < atmos->Nelem;  n++) {
    elemID = atmos->elements[n].ID;
    result &= xdr_counted_string(&xdrs, &elemID);
    result &= xdr_double(&xdrs, &atmos->elements[n].weight);
    result &= xdr_double(&xdrs, &atmos->elements[n].abund);
  }
  /* --- In case magnetic fields were present --       -------------- */

  if (atmos->Stokes) {
    result &= xdr_bool(&xdrs, &atmos->Stokes);

    result &= xdr_vector(&xdrs, (char *) atmos->B, Nspace,
			 sizeof(double), (xdrproc_t) xdr_double);
    result &= xdr_vector(&xdrs, (char *) atmos->gamma_B, Nspace,
			 sizeof(double), (xdrproc_t) xdr_double);
    result &= xdr_vector(&xdrs, (char *) atmos->chi_B, Nspace,
			 sizeof(double), (xdrproc_t) xdr_double);
  }    

  if (!result) {
    sprintf(messageStr, "Unable to write proper amount to output file %s",
	    input.atmos_output);
    Error(ERROR_LEVEL_1, routineName, messageStr);
  }
  xdr_destroy(&xdrs);
  fclose(fp_out);
}
/* ------- end ---------------------------- writeAtmos.c ------------ */
