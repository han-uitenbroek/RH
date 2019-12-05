/* ------- file: -------------------------- writedamp_xdr.c ---------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Mon Jul 29 10:47:04 2002 --

       --------------------------                      ----------RH-- */

/* --- Routines for writing line broadening velocity vbroad and damping
       parameter adamp for each line to file.
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

bool_t xdr_damping(XDR *xdrs, Atom *atom);


/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern InputData input;
extern char messageStr[];


/* ------- begin -------------------------- writeDamping.c ---------- */

bool_t writeDamping(Atom *atom)
{
  const char routineName[] = "writeDamping";

  char  dampfile[15];
  FILE *fp_damp;
  XDR   xdrs;

  /* --- Write broadening velocity and damping parameter -- --------- */

  if (!strcmp(input.dampingFile, "none")) return FALSE;

  sprintf(dampfile, (atom->ID[1] == ' ') ?
	  "damping.%.1s.out" : "damping.%.2s.out", atom->ID);

  if ((fp_damp = fopen(dampfile, "w" )) == NULL) {
    sprintf(messageStr, "Unable to open output file %s", dampfile);
    Error(ERROR_LEVEL_1, routineName, messageStr);
    return FALSE;
  }
  xdrstdio_create(&xdrs, fp_damp, XDR_ENCODE);

  if (!xdr_damping(&xdrs, atom)) {
    sprintf(messageStr, "Unable to write to output file %s", dampfile);
    Error(WARNING, routineName, messageStr);
  }

  xdr_destroy(&xdrs);
  fclose(fp_damp);
  return TRUE;
}
/* ------- end ---------------------------- writeDamping.c -------- */

/* ------- begin -------------------------- readDamping.c --------- */

bool_t readDamping(Atom *atom)
{
  const char routineName[] = "readDamping";

  char  dampfile[15];
  FILE *fp_damp;
  XDR   xdrs;

  if (!strcmp(input.dampingFile, "none")) return FALSE;

  sprintf(dampfile, (atom->ID[1] == ' ') ?
	  "damping.%.1s.out" : "damping.%.2s.out", atom->ID);

  /* --- Read collision rates --                       -------------- */

  if ((fp_damp = fopen(dampfile, "r" )) == NULL) {
    sprintf(messageStr, "Unable to open input file %s", dampfile);
    Error(ERROR_LEVEL_1, routineName, messageStr);
    return FALSE;
  }
  xdrstdio_create(&xdrs, fp_damp, XDR_DECODE);

  if (!xdr_damping(&xdrs, atom)) {
    sprintf(messageStr, "Unable to read from input file %s", dampfile);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }
  xdr_destroy(&xdrs);
  fclose(fp_damp);
  return TRUE;
}
/* ------- end ---------------------------- readDamping.c ----------- */

/* ------- begin -------------------------- xdr_damping.c ----------- */

bool_t xdr_damping(XDR *xdrs, Atom *atom)
{
  register int kr, k;

  bool_t  result = TRUE;
  double *adamp;

  result += xdr_vector(xdrs, (char *) atom->vbroad, atmos.Nspace,
		       sizeof(double), (xdrproc_t) xdr_double);

  /* --- We only export the damping parameters and don't do anything
         on import since we do not store adamp in any structure -- -- */

  if (xdrs->x_op == XDR_ENCODE) {
    adamp = (double *) malloc(atmos.Nspace * sizeof(double));

    for (kr = 0;  kr < atom->Nline;  kr++) {
      if (atom->line[kr].Voigt)
	Damping(&atom->line[kr], adamp);
      else
	for (k = 0;  k < atmos.Nspace;  k++) adamp[k] = 0.0;

      result += xdr_vector(xdrs, (char *) adamp, atmos.Nspace,
			 sizeof(double), (xdrproc_t) xdr_double);
    }
    free(adamp);
  }

  return result;
}
/* ------- end ---------------------------- xdr_damping.c ----------- */
