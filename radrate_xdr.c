/* ------- file: -------------------------- radrate_xdr.c -----------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Tue Sep 17 14:19:48 2002 --

       --------------------------                      ----------RH-- */

/* --- Routines for writing radiative rates to file.
       XDR (external data representation) version. --  -------------- */

 
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "error.h"
#include "inputs.h"
#include "xdr.h"

/* --- Function prototypes --                          -------------- */

bool_t xdr_radrate(XDR *xdrs, Atom *atom);


/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern InputData input;
extern char messageStr[];


/* ------- begin -------------------------- writeRadRate.c - -------- */

bool_t writeRadRate(Atom *atom)
{
  const char routineName[] = "writeRadRate";

  char  ratesfile[15];
  FILE *fp;
  XDR   xdrs;

  /* --- Write radiatve rates for transitions treated in detail.-- -- */

  if (!strcmp(input.radrateFile, "none")) return FALSE;

  sprintf(ratesfile, (atom->ID[1] == ' ') ?
	  "radrate.%.1s.out" : "radrate.%.2s.out", atom->ID);

  if ((fp = fopen(ratesfile, "w" )) == NULL) {
    sprintf(messageStr, "Unable to open output file %s", ratesfile);
    Error(ERROR_LEVEL_1, routineName, messageStr);
    return FALSE;
  }
  xdrstdio_create(&xdrs, fp, XDR_ENCODE);

  if (!xdr_radrate(&xdrs, atom)) {
    sprintf(messageStr, "Unable to write to output file %s", ratesfile);
    Error(WARNING, routineName, messageStr);
  }

  xdr_destroy(&xdrs);
  fclose(fp);
  return TRUE;
}
/* ------- end ---------------------------- writeRadRate.c ---------- */

/* ------- begin -------------------------- readRadRate.c ----------- */

bool_t readRadRate(Atom *atom)
{
  const char routineName[] = "readRadRate";

  char  ratesfile[15];
  FILE *fp;
  XDR   xdrs;

  if (!strcmp(input.radrateFile, "none")) return FALSE;

  sprintf(ratesfile, (atom->ID[1] == ' ') ?
	  "radrate.%.1s.out" : "radrate.%.2s.out", atom->ID);

  /* --- Read radiative rates for transitions treated in detail. -- - */

  if ((fp = fopen(ratesfile, "r" )) == NULL) {
    sprintf(messageStr, "Unable to open input file %s", ratesfile);
    Error(ERROR_LEVEL_2, routineName, messageStr);
    return FALSE;
  }
  xdrstdio_create(&xdrs, fp, XDR_DECODE);

  if (!xdr_radrate(&xdrs, atom)) {
    sprintf(messageStr, "Unable to read from input file %s", ratesfile);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }
  xdr_destroy(&xdrs);
  fclose(fp);
  return TRUE;
}
/* ------- end ---------------------------- readRadRate.c ----------- */

/* ------- begin -------------------------- xdr_radrate.c ----------- */

bool_t xdr_radrate(XDR *xdrs, Atom *atom)
{
  register int kr;

  bool_t result = TRUE;
  AtomicLine *line;
  AtomicContinuum *continuum;

  for (kr = 0;  kr < atom->Nline;  kr++) {
    line = &atom->line[kr];
    result &= xdr_vector(xdrs, (char *) line->Rij, atmos.Nspace,
			 sizeof(double), (xdrproc_t) xdr_double);
    result &= xdr_vector(xdrs, (char *) line->Rji, atmos.Nspace,
			 sizeof(double), (xdrproc_t) xdr_double);
  }
  for (kr = 0;  kr < atom->Ncont;  kr++) {
    continuum = &atom->continuum[kr];
    result &= xdr_vector(xdrs, (char *) continuum->Rij, atmos.Nspace,
			 sizeof(double), (xdrproc_t) xdr_double);
    result &= xdr_vector(xdrs, (char *) continuum->Rji, atmos.Nspace,
			 sizeof(double), (xdrproc_t) xdr_double);
  }
  return result;
}
/* ------- end ---------------------------- xdr_radrate.c ----------- */
