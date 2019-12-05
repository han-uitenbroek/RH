/* ------- file: -------------------------- bessel.c ----------------

       Version:       rh1.0, tools
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Oct 8, 1996

       --------------------------                      ----------RH-- */

/* --- Modified Bessel functions I0, I1, K0, K1.

  See: Numerical Recipes (1st edition), Press et al., p. 176-179
       --                                              -------------- */
 
#include <math.h>

/* --- Function prototypes --                          -------------- */

double  Bessel_I0( double x );
double  Bessel_I1( double x );

/* --- Global variables --                             -------------- */


/* ------- begin -------------------------- Bessel_I0.c ------------- */

double Bessel_I0( double x )
{
  double  ax, ans, y;
  
  if ((ax = fabs(x)) < 3.75) {
    y   = (x * x) / (3.75 * 3.75);
    ans = 1.0 + y*(3.5156229 + y*(3.0899424 + y*(1.2067492 +
	  y*(0.2659732 + y*(0.360768E-1 + y*0.45813E-2)))));
  } else {
    y = 3.75 / ax;
    ans = (exp(ax) / sqrt(ax)) * (0.39894228 + y*(0.1328592E-1 +
	  y*(0.225319E-2 + y*(-0.157565E-2 + y*(0.916281E-2 + 
          y*(-0.2057706E-1 + y*(0.2635537E-1 + y*(-0.1647633E-1 +
          y*0.392377E-2))))))));
  }
  return ans;
}
/* ------- end ---------------------------- Bessel_Io.c ------------- */

/* ------- begin -------------------------- Bessel_I1.c ------------- */

double Bessel_I1( double x )
{
  double  ax, ans, y;
  
  if ((ax = fabs(x)) < 3.75) {
    y   = (x * x) / (3.75 * 3.75);
    ans = ax * (0.5 + y*(0.87890594 + y*(0.51498869 + y*(0.15084934 +
          y*(0.2658733E-1 + y*(0.301532E-2 + y*0.32411E-3))))));
  } else {
    y = 3.75 / ax;
    ans = 0.2282967E-1 + y*(-0.2895312E-1 + y*(0.1787654E-1 - y*0.420059E-2));
    ans = 0.39894228 + y*(-0.3988024E-1 + y*(-0.362018E-2 +
          y*(0.163801E-2 + y*(-0.1031555E-1 + y*ans))));
    ans *= (exp(ax) / sqrt(ax));
  }
  return (x < 0.0) ? -ans : ans;
}
/* ------- end ---------------------------- Bessel_I1.c ------------- */

/* ------- begin -------------------------- Bessel_K0.c ------------- */

double Bessel_K0( double x )
{
  double y, ans;
  
  if (x <= 2.0) {
    y = x * x / 4.0;
    ans = (-log(x / 2.0) * Bessel_I0(x)) + (-0.57721566 + y*(0.42278420 +
	  y*(0.23069756 + y*(0.3488590E-1 + y*(0.262698E-2 +
          y*(0.10750E-3 + y*0.74E-5))))));
  } else {
    y = 2.0 / x;
    ans = (exp(-x) / sqrt(x)) * (1.25331414 + y*(-0.7832358E-1 +
          y*(0.2189568E-1 + y*(-0.1062446E-1 + y*(0.587872E-2  +
          y*(-0.251540E-2 + y*0.53208E-3))))));
  }
  return ans;
}
/* ------- end ---------------------------- Bessel_K0.c ------------- */

/* ------- begin -------------------------- Bessel_K1.c ------------- */

double Bessel_K1( double x )
{
  double y, ans;

  if (x <= 2.0) {
    y = x * x / 4.0;
    ans = (log(x/2.0) * Bessel_I1(x)) + (1.0/x) * (1.0 + y*(0.15443144 +
          y*(-0.67278579 + y*(-0.18156897 + y*(-0.1919402E-1 +
          y*(-0.110404E-2 + y*(-0.4686E-4)))))));
  } else {
    y = 2.0 / x;
    ans = (exp(-x) / sqrt(x)) * (1.25331414 + y*(0.23498619 +
          y*(-0.3655620E-1 + y*(0.1504268E-1 + y*(-0.780353E-2 +
          y*(0.325614E-2 + y*(-0.68245E-3)))))));
  }
  return ans;
}
/* ------- end ---------------------------- Bessel_K1.c ------------- */
