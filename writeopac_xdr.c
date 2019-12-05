/* ------- file: -------------------------- writeopac_xdr.c ---------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Fri May 15 03:09:38 2009 --

       --------------------------                      ----------RH-- */

/* --- Writes opacity and emissivity of active set of transitions for
       each wavelength in the spectrum to file. If atmos.moving == TRUE
       only write the quantities for the upward ray to file.

       Write the record structure of the output file to a separate file
       (ASRS_DOT_OUT). The record structure is an integer array of size 
       spectrum.Nspect in a static atmosphere and spectrum.Nspect *
       atmos.Nrays otherwise. Whenever a wavelength has no overlap with
       an active transition the record number for that wavelength is set
       to -1;

       XDR (external data representation) version. --  -------------- */

 
#include <stdlib.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "spectrum.h"
#include "error.h"
#include "inputs.h"
#include "xdr.h"

#define  ASRS_DOT_OUT  "asrs.out"


/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern Spectrum spectrum;
extern InputData input;
extern char messageStr[];


/* ------- begin -------------------------- writeOpacity.c ---------- */

void writeOpacity(void)
{
  const char routineName[] = "writeOpacity";
  register int nspect, mu;

  bool_t  to_obs, initialize, crosscoupling, boundbound, polarized,
    PRD_angle_dep, result = TRUE;
  int     Nspace = atmos.Nspace, record, *as_rn, Nrecord;
  FILE   *fp_out;
  XDR     xdrs;
  ActiveSet *as;

  if (!strcmp(input.opac_output, "none")) return;

  if ((fp_out = fopen(input.opac_output, "w")) == NULL) {
    sprintf(messageStr, "Unable to open output file %s",
	    input.opac_output);
    Error(ERROR_LEVEL_1, routineName, messageStr);
    return;
  }
  xdrstdio_create(&xdrs, fp_out, XDR_ENCODE);

  if (atmos.moving || atmos.Stokes ||
      (atmos.NPRDactive > 0 && input.PRD_angle_dep))
    Nrecord = atmos.Nrays*spectrum.Nspect;
  else
    Nrecord = spectrum.Nspect;
  as_rn = (int *) malloc(Nrecord * sizeof(int));

  record = -1;
  for (nspect = 0;  nspect < spectrum.Nspect;  nspect++) {
    as = &spectrum.as[nspect];

    if (containsActive(as)) {
      alloc_as(nspect, crosscoupling=FALSE);

      /* --- Check whether current active set includes a bound-bound
             and/or polarized transition and/or angledependent PRD
             transition. Otherwise, only angle-independent opacity and
             source functions are needed --            -------------- */ 

      boundbound    = containsBoundBound(as);
      PRD_angle_dep = (containsPRDline(as) && input.PRD_angle_dep);
      polarized     = containsPolarized(as);

      /* --- Case of angle-dependent opacity and source function -- - */

      if (polarized || PRD_angle_dep || (atmos.moving && boundbound)) {
	for (mu = 0;  mu < atmos.Nrays;  mu++) {
	  initialize = (mu == 0);
	  Opacity(nspect, mu, to_obs=TRUE, initialize);
	  result &= xdr_vector(&xdrs, (char *) as->chi, Nspace, 
			      sizeof(double), (xdrproc_t) xdr_double);
	  result &= xdr_vector(&xdrs, (char *) as->eta, Nspace, 
			      sizeof(double), (xdrproc_t) xdr_double);
	  as_rn[nspect*atmos.Nrays + mu] = ++record; 
	}
      } else {
	Opacity(nspect, 0, to_obs=TRUE, initialize=TRUE);
	result &= xdr_vector(&xdrs, (char *) as->chi, Nspace, 
			    sizeof(double), (xdrproc_t) xdr_double);
	result &= xdr_vector(&xdrs, (char *) as->eta, Nspace, 
			    sizeof(double), (xdrproc_t) xdr_double);
	if (atmos.moving || atmos.Stokes ||
	    (atmos.NPRDactive > 0 && input.PRD_angle_dep)) {
	  record++;
	  for (mu = 0;  mu < atmos.Nrays;  mu++)
	    as_rn[nspect*atmos.Nrays + mu] = record;
	} else
	  as_rn[nspect] = ++record;
      }
      free_as(nspect, crosscoupling=FALSE);
    } else {
      if (atmos.moving || atmos.Stokes ||
	  (atmos.NPRDactive > 0 && input.PRD_angle_dep))
	for (mu = 0;  mu < atmos.Nrays;  mu++)
	  as_rn[nspect*atmos.Nrays + mu] = -1;
      else
	as_rn[nspect] = -1;
    }
  }
  xdr_destroy(&xdrs);
  fclose(fp_out);

  if ((fp_out = fopen(ASRS_DOT_OUT, "w")) == NULL) {
    sprintf(messageStr, "Unable to open output file %s", ASRS_DOT_OUT);
    Error(ERROR_LEVEL_1, routineName, messageStr);
    return;
  }
  xdrstdio_create(&xdrs, fp_out, XDR_ENCODE);

  result &= xdr_vector(&xdrs, (char *) as_rn, Nrecord,
		       sizeof(int), (xdrproc_t) xdr_int);
  free(as_rn);
  xdr_destroy(&xdrs);
  fclose(fp_out);
}
/* ------- end ---------------------------- writeOpacity.c ---------- */
