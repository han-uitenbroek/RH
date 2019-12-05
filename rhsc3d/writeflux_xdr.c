/* ------- file: -------------------------- writeflux_xdr.c ---------

       Version:       rh2.0, 3-D Cartesian, short characteristics
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Tue Nov  7 11:27:59 2000 --

       --------------------------                      ----------RH-- */

 
#include <stdlib.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "geometry.h"
#include "spectrum.h"
#include "constant.h"
#include "error.h"
#include "xdr.h"

/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern Geometry geometry;
extern Spectrum spectrum;
extern char messageStr[];


/* ------- begin -------------------------- writeFlux.c ------------- */

bool_t writeFlux(char *fileName)
{
  const char routineName[] = "writeFlux";
  register int k, mu, nspect;

  bool_t  result = TRUE;
  int     Nrays = geometry.Nrays, Nplane = geometry.Nplane;
  double *flux, *wmuz;
  FILE   *fp_flux;
  XDR     xdrs;

  /* --- Write the radiative flux in the z-direction -- ------------ */

  if ((fp_flux = fopen(fileName, "w")) == NULL) {
    sprintf(messageStr, "Unable to open output file %s", fileName);
    Error(ERROR_LEVEL_1, routineName, messageStr);
    return FALSE;
  }
  xdrstdio_create(&xdrs, fp_flux, XDR_ENCODE);

  flux = (double *) malloc(Nplane * sizeof(double));
  wmuz = (double *) malloc(Nrays * sizeof(double));
  for (mu = 0;  mu < Nrays;  mu++)
    wmuz[mu] = geometry.muz[mu] * geometry.wmu[mu];

  for (nspect = 0;  nspect < spectrum.Nspect;  nspect++) {
    for (k = 0;  k < Nplane;  k++) flux[k] = 0.0;
    for (mu = 0;  mu < geometry.Nrays;  mu++) {
      for (k = 0;  k < Nplane;  k++)
	flux[k] += spectrum.I[nspect*Nrays + mu][k] * wmuz[mu];
    }
    for (k = 0;  k < Nplane;  k++) flux[k] *= 2.0 * PI;
    result &= xdr_vector(&xdrs, (char *) flux, Nplane,
			 sizeof(double), (xdrproc_t) xdr_double);
  }
      
  xdr_destroy(&xdrs);
  fclose(fp_flux);
  free(wmuz);  free(flux);

  return result;
}
/* ------- end ---------------------------- writeFlux.c ------------- */
