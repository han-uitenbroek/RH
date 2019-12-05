/* ------- file: -------------------------- maxchange.c -------------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Tue Jan 25 13:40:59 2000 --

       --------------------------                      ----------RH-- */

/* --- Prints maximum change between iterations using
       the Ng structure. --                            -------------- */

 
#include <math.h>

#include "rh.h"
#include "accelerate.h"
#include "error.h"
#include "inputs.h"


/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

extern CommandLine commandline;
extern char messageStr[];


/* ------- begin -------------------------- MaxChange.c ------------- */

double MaxChange(struct Ng *Ngs, char *text, bool_t quiet)
{
  register int k;

  double dmax = 0.0, *old, *new;

  if (Ngs->count < 2) return dmax;

  old = Ngs->previous[(Ngs->count - 2) % (Ngs->Norder + 2)];
  new = Ngs->previous[(Ngs->count - 1) % (Ngs->Norder + 2)];
  for (k = 0;  k < Ngs->N;  k++) {
    if (new[k])
      dmax = MAX(dmax, fabs((new[k] - old[k]) / new[k]));
  }
  if (!quiet)
    fprintf(commandline.logfile, "%s delta = %6.4E", text, dmax);

  return dmax;
}
/* ------- end ---------------------------- MaxChange.c ------------- */
