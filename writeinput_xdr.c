/* ------- file: -------------------------- writeinput_xdr.c --------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Wed Apr  1 09:43:12 2009 --

       --------------------------                      ----------RH-- */

/* --- XDR (external data representation) version. --  -------------- */


#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "spectrum.h"
#include "error.h"
#include "inputs.h"
#include "xdr.h"

#define  INPUT_DOT_OUT  "input.out"


/* --- Function prototypes --                          -------------- */

int is_big_endian(void);


/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern Spectrum spectrum;
extern InputData input;
extern char messageStr[];


/* ------- begin -------------------------- writeInput.c ------------ */

void writeInput(void)
{
  const char routineName[] = "writeInput";

  bool_t  result=TRUE, PRD_angle_dep, XRD, big_endian;
  FILE   *fp_out;
  XDR     xdrs;

  if (!strcmp(INPUT_DOT_OUT, "none")) return;

  if ((fp_out = fopen(INPUT_DOT_OUT, "w")) == NULL) {
    sprintf(messageStr, "Unable to open output file %s",
	    INPUT_DOT_OUT);
    Error(ERROR_LEVEL_1, routineName, messageStr);
    return;
  }
  xdrstdio_create(&xdrs, fp_out, XDR_ENCODE);

  PRD_angle_dep = (input.PRD_angle_dep  &&  atmos.NPRDactive > 0);
  XRD           = (input.XRD  &&  atmos.NPRDactive > 0);

  /* --- Write various input parameters to file --     -------------- */      

  result &= xdr_bool(&xdrs, &input.magneto_optical);
  result &= xdr_bool(&xdrs, &PRD_angle_dep);
  result &= xdr_bool(&xdrs, &XRD);

  result &= xdr_enum(&xdrs, (enum_t *) &input.startJ);
  result &= xdr_enum(&xdrs, (enum_t *) &input.StokesMode);

  result &= xdr_double(&xdrs, &input.metallicity);

  result &= xdr_bool(&xdrs, &input.backgr_pol);

  /* --- Write Endianness of compute architecture so that J can be
         read properly in the analysis --              -------------- */

  big_endian = is_big_endian();
  result &= xdr_bool(&xdrs, &big_endian);


  if (!result) {
    sprintf(messageStr, "Unable to write proper amount to output file %s",
	    INPUT_DOT_OUT);
    Error(ERROR_LEVEL_1, routineName, messageStr);
  }
  xdr_destroy(&xdrs);
  fclose(fp_out);
}
/* ------- end ---------------------------- writeInput.c ------------ */

/* ------- begin -------------------------- is_big_endian.c --------- */

int is_big_endian(void)
{
  unsigned char *b;
  short i[1] = {1};

  /* --- Returns 1 if host machine is big endian, 0 otherwise -- ---- */

  b = (unsigned char *) i;
  return (b[0]) ? 0 : 1;
}
/* ------- end ---------------------------- is_big_endian.c --------- */
