/* ------- file: -------------------------- xdr.h -------------------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Tue Sep 21 13:36:03 2004 --

       --------------------------                      ----------RH-- */

#ifndef __XDR_H__
#define __XDR_H__

/* --- To be included in routines that use the XDR (see ``man xdr'').
       (eXternal Data Representation) format to read/write data.
       --                                              -------------- */


#include <rpc/xdr.h>

#include "atom.h"
#include "atmos.h"

/* --- Associated function prototypes --               -------------- */

void   convertJ(FILE *fp_J, int Nspect, enum xdr_op xdrOperation);
bool_t readB(Atmosphere *atmos);
bool_t writeFlux(char *fileName);
bool_t writeMetals(char *fileName);
bool_t xdr_atom(XDR *xdrs, Atom *atom);
bool_t xdr_counted_string(XDR *xdrs, char **p);
bool_t xdr_populations(XDR *xdrs, char *atmosID, int Nlevel, int Nspace,
		       double *n, double *nstar);
bool_t xdr_BRS(XDR *xdrs);


#endif /* !__XDR_H__ */

/* ------- end ---------------------------- xdr.h ------------------- */
