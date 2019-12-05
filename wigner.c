/* ------- file: -------------------------- wigner.c ----------------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Mon Aug  3 15:52:07 2009   --

       --------------------------                      ----------RH-- */

/* --- Functions for evaluation of the Wigner 3-J, 6-J symbols -- --- */

#include <stdlib.h>
#include <math.h>

#include "rh.h"
#include "atom.h"


#define NFACT 301


/* --- Function prototypes --                          -------------- */

double fact(double xi);
double triangleCoefficient(double a, double b, double c);


/* --- Global variables --                             -------------- */


/* ------- begin -------------------------- w3js.c ------------------ */

double w3js(double j1, double j2, double j3,
	    double m1, double m2, double m3)
{
  register int t;

  int    tmin, tmax;
  double w3js, q1, q2, q3, q4, q5, sumt;

  /* --- See: E. Landi Degl'Innocenti & M. Landolfi, 2004
              "Polarization in Spectral Lines", Kluwer

         Eqs. 2.19 and 2.22

         See also: http://mathworld.wolfram.com/Wigner3j-Symbol.html

         --                                            -------------- */

   w3js = 0.0;

  /* --- Check whether j1, ..., m3 are whole or half-integers -- ---- */

  if (fmod(j1, 0.5) || fmod(m1, 0.5) ||
      fmod(j2, 0.5) || fmod(m2, 0.5) ||
      fmod(j3, 0.5) || fmod(m3, 0.5)) return w3js;

  /* --- Check whether the m_i are between -j_i and j_i -- ---------- */
							 
  if (fabs(m1) > j1  ||  fabs(m2) > j2  ||  fabs(m3) > j3) return w3js;

  /* --- Check sum rule for z-components of J --       -------------- */

  if (m1 + m2 + m3 != 0.0) return w3js;

  /* --- Check maximum value of j3, given j1,2 --      -------------- */

  if (j3 > j1 + j2  ||  j3 < fabs(j1 - j2)) return w3js;

  q1 = j3 - j2 + m1;
  q2 = j3 - j1 - m2;
  tmin = (int) MAX(0, MAX(-q1, -q2));

  q3 = j1 + j2 - j3;
  q4 = j1 - m1;
  q5 = j2 + m2;
  tmax = (int) MIN(q3, MIN(q4, q5));

  sumt = 0.0;
  for (t = tmin;  t <= tmax;  t++) {
    sumt += ((t % 2) ? -1.0 : 1.0) /
      (fact(t) * fact(q3 - t) * fact(q4 - t) * fact(q5 - t) *
       fact(q1 + t) * fact(q2 + t));
  }
  
  w3js = ((((int) (j1 - j2 - m3)) % 2) ? -1.0 : 1.0) *
    sqrt(triangleCoefficient(j1, j2, j3) *
	 fact(j1 + m1) * fact(j1 - m1) *
	 fact(j2 + m2) * fact(j2 - m2) *
	 fact(j3 + m3) * fact(j3 - m3)) * sumt;

  return w3js;
}
/* ------- end ---------------------------- w3js.c ------------------ */

/* ------- begin -------------------------- w6js.c ------------------ */

double w6js(double j1, double j2, double j3,
	    double J1, double J2, double J3)
{
  register int z;

  int    zmin, zmax;
  double w6js, q1, q2, q3, sum1, sum2, sum3, sum4,
    tc1, tc2, tc3, tc4, sumz;

  /* --- See: E. Landi Degl'Innocenti & M. Landolfi, 2004
              "Polarization in Spectral Lines", Kluwer

         Eq. 2.35

         See also: http://mathworld.wolfram.com/Wigner6j-Symbol.html

         --                                            -------------- */

  w6js = 0.0;

  /* --- Check whether j1, ..., m3 are whole or half-integers -- ---- */

  if (fmod(j1, 0.5) || fmod(J1, 0.5) ||
      fmod(j2, 0.5) || fmod(J2, 0.5) ||
      fmod(j3, 0.5) || fmod(J3, 0.5)) return w6js;

  if (j3 > j1 + j2  ||  j3 < fabs(j1 - j2)) return w6js;
  if (J3 > j1 + J2  ||  J3 < fabs(j1 - J2)) return w6js;
  if (J3 > J1 + j2  ||  J3 < fabs(J1 - j2)) return w6js;
  if (j3 > J1 + J2  ||  j3 < fabs(J1 - J2)) return w6js;

  if (fmod(sum1 = j1 + j2 + j3, 1.0)) return w6js;
  if (fmod(sum2 = j1 + J2 + J3, 1.0)) return w6js;
  if (fmod(sum3 = J1 + j2 + J3, 1.0)) return w6js;
  if (fmod(sum4 = J1 + J2 + j3, 1.0)) return w6js;
  
  zmin = MAX(MAX(MAX(sum1, sum2), sum3), sum4);

  q1 = j1 + j2 + J1 + J2;
  q2 = j2 + j3 + J2 + J3;
  q3 = j3 + j1 + J3 + J1;
  zmax = MIN(MIN(q1, q2), q3);

  sumz = 0.0;
  for (z = zmin;  z <= zmax;  z++) {
    sumz += ((z % 2) ? -1.0 : 1.0) * fact(z + 1.0) /
      (fact(z - sum1) * fact(z - sum2) * fact(z - sum3) * fact(z - sum4) *
       fact(q1 - z) * fact(q2 - z) * fact(q3 - z));
  }

  tc1 = triangleCoefficient(j1, j2, j3);
  tc2 = triangleCoefficient(j1, J2, J3);
  tc3 = triangleCoefficient(J1, j2, J3);
  tc4 = triangleCoefficient(J1, J2, j3);

  w6js = sqrt(tc1 * tc2 * tc3 * tc4) * sumz;

  return w6js;
}
/* ------- begin -------------------------- fact.c ------------------ */

double fact(double xi)
{
  /* --- Returns the factorial of xi. On first call the factorials
         for the numbers 0 through 300 are stored in the static
         variable factorial. --                        -------------- */

  register int n;

  static bool_t initialize = TRUE;
  static double *factorial;

  if (initialize) {
    factorial = (double *) malloc(NFACT * sizeof(double));

    factorial[0]  = 1.0;
    for (n = 1;  n < NFACT;  n++) factorial[n] = factorial[n-1] * n;
 
    initialize = FALSE;
  }
  /* --- Test whether xi is an integer, otherwise return 1.0 -- ----- */

  return (fmod(xi, 1.0)) ? factorial[0] : factorial[(int) xi];
}
/* ------- end ---------------------------- fact.c ------------------ */

/* ------- begin -------------------------- triangleCoefficient.c --- */

double triangleCoefficient(double a, double b, double c)
{
  /* --- Returns the square of the triangular coefficient 
         \Delta(a, b, c). Square root is not applied for efficiency
         in calling program.

         See: http://mathworld.wolfram.com/TriangleCoefficient.html

         --                                            -------------- */

  double tricoeff;

  tricoeff = fact(a + b - c) * fact(a - b + c) * fact(-a + b + c) /
    fact(a + b + c + 1.0);

  return tricoeff;
}
/* ------- end ---------------------------- triangleCoefficient.c --- */
