/* ------- file: -------------------------- writemetal_xdr.c --------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Tue Apr 21 15:39:29 2009 --

       --------------------------                      ----------RH-- */

/* --- Writes atomic data of background opacity package to output file.
       XDR (external data representation) version. --  -------------- */

 
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "error.h"
#include "xdr.h"

/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern char messageStr[];


/* ------- begin -------------------------- writeMetals.c ----------- */

bool_t writeMetals(char *fileName)
{
  const char routineName[] = "writeMetals";
  register int n;

  int    result = TRUE, Nmetal;
  Atom *atom;
  FILE  *fp_out;
  XDR    xdrs;

  if (!strcmp(fileName, "none")) return result;

  if ((fp_out = fopen(fileName, "w")) == NULL) {
    sprintf(messageStr, "Unable to open output file %s", fileName);
    Error(ERROR_LEVEL_1, routineName, messageStr);
    return FALSE;
  }
  xdrstdio_create(&xdrs, fp_out, XDR_ENCODE);

  xdr_atom(&xdrs, atmos.H);

  /* --- Only write active atoms, write Hydrogen seperately -- ------ */

  Nmetal = atmos.Natom - 1 - atmos.Nactiveatom;
  result &= xdr_int(&xdrs, &Nmetal);
  for (n = 1;  n < atmos.Natom;  n++) {
    atom = &atmos.atoms[n];
    if (!atom->active) 
      result &= xdr_atom(&xdrs, atom);
  }

  xdr_destroy(&xdrs);
  fclose(fp_out);

  return result;
}
/* ------- end ---------------------------- writeMetals.c ----------- */
