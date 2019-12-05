/* ------- file: -------------------------- brs_xdr.c ---------------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Thu Jan  3 14:35:46 2008 --

       --------------------------                      ----------RH-- */

/* --- Routines to write and read background record structure to file.

       XDR (external data representation) version. --  -------------- */


#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "spectrum.h"
#include "error.h"
#include "xdr.h"

#define  BRS_DOT_OUT  "brs.out"


/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern Spectrum spectrum;
extern char messageStr[];


/* ------- begin -------------------------- writeBRS.c -------------- */

void writeBRS(void)
{
  const char routineName[] = "writeBRS";

  FILE *fp_out;
  XDR   xdrs;


  if ((fp_out = fopen(BRS_DOT_OUT, "w")) == NULL) {
    sprintf(messageStr, "Unable to open output file %s", BRS_DOT_OUT);
    Error(ERROR_LEVEL_1, routineName, messageStr);
    return;
  }
  xdrstdio_create(&xdrs, fp_out, XDR_ENCODE);

  if (!xdr_BRS(&xdrs)) {
    sprintf(messageStr, "Unable to write to output file %s", BRS_DOT_OUT);
    Error(ERROR_LEVEL_1, routineName, messageStr);
  }

  xdr_destroy(&xdrs);
  fclose(fp_out);
}
/* ------- end ---------------------------- writeBRS.c -------------- */

/* ------- begin -------------------------- readBRS.c --------------- */

void readBRS(void)
{
  const char routineName[] = "readBRS";

  FILE *fp_in;
  XDR   xdrs;

  if ((fp_in = fopen(BRS_DOT_OUT, "r")) == NULL) {
    sprintf(messageStr, "Unable to open input file %s", BRS_DOT_OUT);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }
  xdrstdio_create(&xdrs, fp_in, XDR_DECODE);

  if (!xdr_BRS(&xdrs)) {
    sprintf(messageStr, "Unable to read from input file %s", BRS_DOT_OUT);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }

  xdr_destroy(&xdrs);
  fclose(fp_in);
}
/* ------- end ---------------------------- readBRS.c --------------- */

/* ------- begin -------------------------- xdr_BRS.c --------------- */

bool_t xdr_BRS(XDR *xdrs)
{
  const char routineName[] = "xdr_BRS";
  register int nspect;

  char  *atmosID;
  bool_t result = TRUE, *hasline, *ispolarized;
  int    Nspace, Nsp, Nrecno;

  hasline     =
    (bool_t *) malloc(spectrum.Nspect * sizeof(bool_t));
  ispolarized =
    (bool_t *) malloc(spectrum.Nspect * sizeof(bool_t));

  if (atmos.moving || atmos.Stokes)
    Nrecno = 2 * spectrum.Nspect * atmos.Nrays;
  else
    Nrecno = spectrum.Nspect;

  if (xdrs->x_op == XDR_ENCODE) {
    atmosID = atmos.ID;
    result &= xdr_counted_string(xdrs, &atmosID);
    result &= xdr_int(xdrs, &atmos.Nspace);
    result &= xdr_int(xdrs, &spectrum.Nspect);

    /* --- Flags for presence of background line --    -------------- */
 
    for (nspect = 0;  nspect < spectrum.Nspect;  nspect++)
      hasline[nspect] = atmos.backgrflags[nspect].hasline;
    result &= xdr_vector(xdrs, (char *) hasline, spectrum.Nspect,
			 sizeof(bool_t), (xdrproc_t) xdr_bool);

    /* --- Flags for presence of polarized background -- ------------ */
 
    for (nspect = 0;  nspect < spectrum.Nspect;  nspect++)
      ispolarized[nspect] = atmos.backgrflags[nspect].ispolarized;
    result &= xdr_vector(xdrs, (char *) ispolarized, spectrum.Nspect,
			 sizeof(bool_t), (xdrproc_t) xdr_bool);
  } else {
    atmos.backgrflags = (flags *) malloc(spectrum.Nspect * sizeof(flags));
    atmos.backgrrecno = (int *) malloc(Nrecno * sizeof(int));

    result &= xdr_counted_string(xdrs, &atmosID);
    if (!strstr(atmosID, atmos.ID)) {
      sprintf(messageStr,
	      "Input was derived from different atmosphere (%s) than current",
              atmosID);
      Error(WARNING, routineName, messageStr);
    }
    free(atmosID);

    result &= xdr_int(xdrs, &Nspace);
    result &= xdr_int(xdrs, &Nsp);
    if (Nspace != atmos.Nspace || Nsp != spectrum.Nspect) {
      free(hasline);
      return FALSE;
    }
    /* --- Flags for presence of background line --    -------------- */
 
    result &= xdr_vector(xdrs, (char *) hasline, spectrum.Nspect,
			 sizeof(bool_t), (xdrproc_t) xdr_bool);
    for (nspect = 0;  nspect < spectrum.Nspect;  nspect++)
      atmos.backgrflags[nspect].hasline = hasline[nspect];

    /* --- Flags for presence of polarized background -- ------------ */
 
    result &= xdr_vector(xdrs, (char *) ispolarized, spectrum.Nspect,
			 sizeof(bool_t), (xdrproc_t) xdr_bool);
    for (nspect = 0;  nspect < spectrum.Nspect;  nspect++)
      atmos.backgrflags[nspect].ispolarized = ispolarized[nspect];
  }

  result &= xdr_vector(xdrs, (char *) atmos.backgrrecno, Nrecno,
		       sizeof(int), (xdrproc_t) xdr_int);

  free(hasline);
  free(ispolarized);

  return result;
}
/* ------- end ---------------------------- xdr_BRS.c --------------- */
