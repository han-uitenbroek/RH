/* ------- file: -------------------------- writecoll_xdr.c ---------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Fri Jul 26 13:38:58 2002 --

       --------------------------                      ----------RH-- */

/* --- Routines for writing collisional rates to file.
       XDR (external data representation) version. --  -------------- */

 
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "error.h"
#include "inputs.h"
#include "xdr.h"

/* --- Function prototypes --                          -------------- */

bool_t xdr_collrate(XDR *xdrs, Atom *atom);


/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern InputData input;
extern char messageStr[];


/* ------- begin -------------------------- writeCollisionRate.c ---- */

bool_t writeCollisionRate(Atom *atom)
{
  const char routineName[] = "writeCollisionRate";

  char  ratesfile[16];
  FILE *fp_rate;
  XDR   xdrs;

  /* --- Write collision rates --                      -------------- */

  if (!strcmp(input.collrateFile, "none")) return FALSE;

  sprintf(ratesfile, (atom->ID[1] == ' ') ?
	  "collrate.%.1s.out" : "collrate.%.2s.out", atom->ID);

  if ((fp_rate = fopen(ratesfile, "w" )) == NULL) {
    sprintf(messageStr, "Unable to open output file %s", ratesfile);
    Error(ERROR_LEVEL_1, routineName, messageStr);
    return FALSE;
  }
  xdrstdio_create(&xdrs, fp_rate, XDR_ENCODE);

  if (!xdr_collrate(&xdrs, atom)) {
    sprintf(messageStr, "Unable to write to output file %s", ratesfile);
    Error(WARNING, routineName, messageStr);
  }

  xdr_destroy(&xdrs);
  fclose(fp_rate);
  return TRUE;
}
/* ------- end ---------------------------- writeCollisionRate.c ---- */

/* ------- begin -------------------------- readCollisionRate.c ----- */

bool_t readCollisionRate(Atom *atom)
{
  const char routineName[] = "readCollisionRate";

  char  ratesfile[16];
  FILE *fp_rate;
  XDR   xdrs;

  if (!strcmp(input.collrateFile, "none")) return FALSE;

  sprintf(ratesfile, (atom->ID[1] == ' ') ?
	  "collrate.%.1s.out" : "collrate.%.2s.out", atom->ID);

  /* --- Read collision rates --                       -------------- */

  if ((fp_rate = fopen(ratesfile, "r" )) == NULL) {
    sprintf(messageStr, "Unable to open input file %s", ratesfile);
    Error(ERROR_LEVEL_1, routineName, messageStr);
    return FALSE;
  }
  xdrstdio_create(&xdrs, fp_rate, XDR_DECODE);

  if (!xdr_collrate(&xdrs, atom)) {
    sprintf(messageStr, "Unable to read from input file %s", ratesfile);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }
  xdr_destroy(&xdrs);
  fclose(fp_rate);
  return TRUE;
}
/* ------- end ---------------------------- readCollisionRate.c ----- */

/* ------- begin -------------------------- xdr_collrate.c ---------- */

bool_t xdr_collrate(XDR *xdrs, Atom *atom)
{
  register int ij;

  bool_t result = TRUE;

  for (ij = 0;  ij < SQ(atom->Nlevel);  ij++) {
    result &= xdr_vector(xdrs, (char *) atom->C[ij], atmos.Nspace,
			 sizeof(double), (xdrproc_t) xdr_double);
  }
  return result;
}
/* ------- end ---------------------------- xdr_collrate.c ---------- */
