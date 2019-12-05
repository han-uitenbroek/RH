/* ------- file: -------------------------- writeatom_xdr.c ---------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Tue Apr 21 16:01:09 2009 --

       --------------------------                      ----------RH-- */

/* --- Writes atomic data to output file.
       XDR (external data representation) version. --  -------------- */

 
#include <stdlib.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "spectrum.h"
#include "inputs.h"
#include "constant.h"
#include "error.h"
#include "xdr.h"


/* --- Function prototypes --                          -------------- */

bool_t xdr_adamp(XDR *xdrs, struct Atom *atom);


/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern Spectrum spectrum;
extern InputData input;
extern char messageStr[];


/* ------- begin -------------------------- writeAtom.c ------------- */

void writeAtom(struct Atom *atom)
{
  const  char routineName[] = "writeAtom";

  char  atomoutfile[12];
  FILE  *fp_out, *fp_adamp;
  XDR    xdrs;

  sprintf(atomoutfile, (atom->ID[1] == ' ') ?
	  "atom.%.1s.out" : "atom.%.2s.out", atom->ID);

  if ((fp_out = fopen(atomoutfile, "w")) == NULL) {
    sprintf(messageStr, "Unable to open output file %s",
	    atomoutfile);
    Error(ERROR_LEVEL_1, routineName, messageStr);
    return;
  }
  xdrstdio_create(&xdrs, fp_out, XDR_ENCODE);

  if (!xdr_atom(&xdrs, atom)) {
    sprintf(messageStr, "Unable to write to output file %s",
	    atomoutfile);
    Error(ERROR_LEVEL_1, routineName, messageStr);
  }
  xdr_destroy(&xdrs);
  fclose(fp_out);
}
/* ------- end ---------------------------- writeAtom.c ------------- */

/* ------- begin -------------------------- xdr_atom.c -------------- */

bool_t xdr_atom(XDR *xdrs, struct Atom *atom)
{
  register int i, kr, kf;

  bool_t result = TRUE, shape;
  enum type rt_type;
  double lambda_air;
  AtomicLine *line;
  AtomicContinuum *continuum;
  FixedTransition *ft;

  /* --- Does the actual writing --                    -------------- */

  result &= xdr_bool(xdrs, &atom->active);

  result &= xdr_int(xdrs, &atom->Nlevel);
  result &= xdr_int(xdrs, &atom->Nline);
  result &= xdr_int(xdrs, &atom->Ncont);
  result &= xdr_int(xdrs, &atom->Nfixed);
  
  result &= xdr_double(xdrs, &atom->abundance);
  result &= xdr_double(xdrs, &atom->weight);
  
  for (i = 0;  i < atom->Nlevel;  i++)
    result &= xdr_counted_string(xdrs, &atom->label[i]);

  result &= xdr_vector(xdrs, (char *) atom->g, atom->Nlevel,
		      sizeof(double), (xdrproc_t) xdr_double);
  result &= xdr_vector(xdrs, (char *) atom->E, atom->Nlevel,
		      sizeof(double), (xdrproc_t) xdr_double);
  result &= xdr_vector(xdrs, (char *) atom->stage, atom->Nlevel,
		      sizeof(int), (xdrproc_t) xdr_int);

  /* --- Write the bound-bound and bound-free transitions -- -------- */

  rt_type = ATOMIC_LINE;
  for (kr = 0;  kr < atom->Nline;  kr++) {
    line = atom->line + kr;

    result &= xdr_enum(xdrs, (enum_t *) &rt_type);
    result &= xdr_int(xdrs, &line->i);
    result &= xdr_int(xdrs, &line->j);
    result &= xdr_int(xdrs, &line->Nlambda);
    result &= xdr_int(xdrs, &line->Nblue);

    if (spectrum.vacuum_to_air) {
      vacuum_to_air(1, &line->lambda0, &lambda_air);
      result &= xdr_double(xdrs, &lambda_air);
    } else
      result &= xdr_double(xdrs, &line->lambda0);

    if (line->Voigt) {
      if (line->PRD) shape = 2;  else  shape = 1;
    } else
      shape = 0;
    result &= xdr_int(xdrs, &shape);
    result &= xdr_double(xdrs, &line->Aji);
  }
  rt_type = ATOMIC_CONTINUUM;
  for (kr = 0;  kr < atom->Ncont;  kr++) {
    continuum = atom->continuum + kr;

    result &= xdr_enum(xdrs, (enum_t *) &rt_type);
    result &= xdr_int(xdrs, &continuum->i);
    result &= xdr_int(xdrs, &continuum->j);
    result &= xdr_int(xdrs, &continuum->Nlambda);
    result &= xdr_int(xdrs, &continuum->Nblue);

    if (spectrum.vacuum_to_air) {
      vacuum_to_air(1, &continuum->lambda0, &lambda_air);
      result &= xdr_double(xdrs, &lambda_air);
    } else
      result &= xdr_double(xdrs, &continuum->lambda0);

    if (continuum->hydrogenic)
      shape = 3;
    else
      shape = 4;
    result &= xdr_int(xdrs, &shape);
    result &= xdr_double(xdrs, &continuum->alpha0);
  }
  /* --- Write wavelength dependent cross section for bound-free -- - */

  for (kr = 0;  kr < atom->Ncont;  kr++) {
    continuum = atom->continuum + kr;

    if (continuum->hydrogenic) {
      result &= xdr_double(xdrs, continuum->lambda);
    } else {
      result &= xdr_vector(xdrs, (char *) continuum->lambda,
			   continuum->Nlambda,
			   sizeof(double), (xdrproc_t) xdr_double);
      result &= xdr_vector(xdrs, (char *) continuum->alpha,
			   continuum->Nlambda,
			   sizeof(double), (xdrproc_t) xdr_double);
    }
  }
  /* --- Write the fixed transtitions --               -------------- */

  for (kf = 0;  kf < atom->Nfixed;  kf++) {
    ft = atom->ft + kf;

    result &= xdr_enum(xdrs, (enum_t *) &ft->type);
    result &= xdr_enum(xdrs, (enum_t *) &ft->option);
    result &= xdr_int(xdrs, &ft->i);
    result &= xdr_int(xdrs, &ft->j);

    if (spectrum.vacuum_to_air) {
      vacuum_to_air(1, &ft->lambda0, &lambda_air);
      result &= xdr_double(xdrs, &lambda_air);
    } else
      result &= xdr_double(xdrs, &ft->lambda0);

    result &= xdr_double(xdrs, &ft->strength);
    result &= xdr_double(xdrs, &ft->Trad);
  }

  return result;
}
/* ------- end ---------------------------- xdr_atom.c -------------- */

/* ------- begin -------------------------- xdr_adamp.c ------------- */

bool_t xdr_adamp(XDR *xdrs, Atom *atom)
{
  register int kr, k;

  int     result = TRUE;
  double *adamp;
  AtomicLine *line;

  adamp = (double *) malloc(atmos.Nspace * sizeof(double));

  for (kr = 0;  kr < atom->Nline;  kr++) {
    line = atom->line + kr;
    if (line->Voigt)
      Damping(line, adamp);
    else
      for (k = 0;  k < atmos.Nspace;  k++) adamp[k] = 0.0;

    result &= xdr_vector(xdrs, (char *) adamp, atmos.Nspace,
			 sizeof(double), (xdrproc_t) xdr_double);
  }
  free(adamp);
  return result;
}
/* ------- end ---------------------------- xdr_adamp.c ------------- */
