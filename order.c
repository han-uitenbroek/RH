/* ------- file: -------------------------- order.c -----------------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Thu Feb  4 22:09:06 1999 --

       --------------------------                      ----------RH-- */

/* --- Functions for ordering arrays of values. See ``man qsort''
       for examples. --                                -------------- */

 
#include "rh.h"

/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */


/* ------- begin -------------------------- qsascend.c -------------- */

int qsascend(const void *v1, const void *v2 )
{
  double f1, f2;

  f1 = *(double*) v1;  f2 = *(double *) v2;
  if (f1 < f2)
    return -1;
  else if (f1 > f2)
    return 1;
  else
    return 0;
}
/* ------- end ---------------------------- qsascend.c -------------- */

/* ------- begin -------------------------- qsdescend.c ------------- */

int qsdescend(const void *v1, const void *v2 )
{
  double f1, f2;

  f1 = *(double*) v1;  f2 = *(double *) v2;
  if (f1 > f2)
    return -1;
  else if (f1 < f2)
    return 1;
  else
    return 0;
}
/* ------- end ---------------------------- qsdescend.c ------------- */
