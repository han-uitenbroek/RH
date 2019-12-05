/* ------- file: -------------------------- gammafunc.c -------------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Mon Sep 12 22:15:56 2005 --

       --------------------------                      ----------RH-- */

/* --- Logarithm of the Gamma function --              -------------- */

#include <math.h>

#include "rh.h"

/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */


/* ------- begin -------------------------- gammln.c ---------------- */

double gammln(double xx)
{
  register int j;

  static double cof[6] = { 76.18009172947146,
			  -86.50532032941677,
			   24.01409824083091,
			   -1.231739572450155,
			    0.1208650973866179e-2,
			   -0.5395239384953e-5};

  double x, y, tmp, ser;

  x = xx;
  y = x;
  tmp = x + 5.5;
  tmp -= (x + 0.5) * log(tmp);
  ser = 1.000000000190015;

  for (j = 0;  j < 6;  j++)
    ser += cof[j] / ++y;

  return -tmp + log(2.5066282746310005 * ser / x);
}
/* ------- end ---------------------------- gammln.c ---------------- */
