/* ------- file: -------------------------- writemolec_xdr.c --------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Wed Sep 22 13:51:03 2004 --

       --------------------------                      ----------RH-- */

/* --- Write and read molecular data of opacity package to output file.
       XDR (external data representation) version. --  -------------- */

 
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "error.h"
#include "xdr.h"


/* --- Function prototypes --                          -------------- */

bool_t xdr_molecules(XDR *xdrs, Molecule *molecules);


/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern char messageStr[];


/* ------- begin -------------------------- writeMolecules.c -------- */

void writeMolecules(char *fileName)
{
  const char routineName[] = "writeMolecules";

  FILE  *fp_molecules;
  XDR    xdrs;

  if (!strcmp(fileName, "none")) return;

  if ((fp_molecules = fopen(fileName, "w")) == NULL) {
    sprintf(messageStr, "Unable to open output file %s", fileName);
    Error(ERROR_LEVEL_1, "writeMolecules", messageStr);
  }
  xdrstdio_create(&xdrs, fp_molecules, XDR_ENCODE);

  if (!xdr_molecules(&xdrs, atmos.molecules)) {
    sprintf(messageStr, "Unable to write to output file %s", fileName);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }

  xdr_destroy(&xdrs);
  fclose(fp_molecules);
}
/* ------- end ---------------------------- writeMolecules.c -------- */

/* ------- begin -------------------------- readMolecules.c --------- */

void readMolecules(char *fileName)
{
  const char routineName[] = "readMolecules";

  FILE  *fp_molecules;
  XDR    xdrs;

  if (!strcmp(fileName, "none")) return;

  if ((fp_molecules = fopen(fileName, "r")) == NULL) {
    sprintf(messageStr, "Unable to open intput file %s", fileName);
    Error(ERROR_LEVEL_1, "readMolecules", messageStr);
  }
  xdrstdio_create(&xdrs, fp_molecules, XDR_DECODE);

  if (!xdr_molecules(&xdrs, atmos.molecules)) {
    sprintf(messageStr, "Unable to read from input file %s", fileName);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }

  xdr_destroy(&xdrs);
  fclose(fp_molecules);
}
/* ------- end ---------------------------- readMolecules.c --------- */

/* ------- begin -------------------------- xdr_molecules.c --------- */

bool_t xdr_molecules(XDR *xdrs, Molecule *molecules)
{
  register int n, kr;
  bool_t result = TRUE;

  char    *molecID;
  double **E;
  Molecule *molecule;
  MolecularLine *mrt;

  result &= xdr_int(xdrs, &atmos.Nmolecule);

  for (n = 0;  n < atmos.Nmolecule;  n++) {
    molecule = &molecules[n];

    molecID = molecule->ID;
    result &= xdr_counted_string(xdrs, &molecID);
    result &= xdr_int(xdrs, &molecule->Nv);
    result &= xdr_int(xdrs, &molecule->NJ);

    result &= xdr_double(xdrs, &molecule->Ediss);
    result &= xdr_double(xdrs, &molecule->Tmin);
    result &= xdr_double(xdrs, &molecule->Tmax);

    result &= xdr_vector(xdrs, (char *) molecule->n, atmos.Nspace,
			 sizeof(double), (xdrproc_t) xdr_double);

    if (molecule->Nv > 0  &&  molecule->NJ > 0) {
      E = matrix_double(molecule->Nv, molecule->NJ);

      if (xdrs->x_op == XDR_ENCODE) {
	for (kr = 0;  kr < molecule->Nrt;  kr++) {
	  mrt = molecule->mrt + kr;

	  E[mrt->vi][(int) (mrt->gi - 1)/2] = mrt->Ei;
	  E[mrt->vj][(int) (mrt->gj - 1)/2] = mrt->Ej;
	}
      }
      result &= xdr_vector(xdrs, (char *) E[0], molecule->Nv * molecule->NJ,
			   sizeof(double), (xdrproc_t) xdr_double);
      freeMatrix((void **) E);
    }
  }
  result &= xdr_vector(xdrs, (char *) atmos.nHmin, atmos.Nspace,
		       sizeof(double), (xdrproc_t) xdr_double);

  return result;
}
/* ------- end ---------------------------- xdr_molecules.c --------- */

