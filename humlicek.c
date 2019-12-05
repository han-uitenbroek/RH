/* ------- file: -------------------------- humlicek.c --------------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Thu Dec 13 17:17:04 2001 --

       --------------------------                      ----------RH-- */

#include "rh.h"
#include "complex.h"

/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

/* -- Voigt function subroutines in different regions.
      Relative accuracy 1.0E-04. Also calculates Faraday-Voigt
      function needed in Stokes radiative transfer.

      C version.

 See: Humlicek 1982, JQSRT 27, p. 437
      --                                               -------------- */

/* ------- begin -------------------------- Humlicek1.c ------------- */

complex Humlicek1(complex z)
{
  /* --- Approximation in region I --                  -------------- */

  complex z1, z2;

  z1 = cmplx_sclr(0.5641896, z);
  z2 = cmplx_addr(cmplx_mult(z, z), 0.5);

  return cmplx_div(z1, z2);
}
/* ------- end ---------------------------- Humlicek1.c ------------- */

/* ------- begin -------------------------- Humlicek2.c ------------- */

complex Humlicek2(complex z)
{
  /* --- Approximation in region II --                 -------------- */

  complex z1, z2, u = cmplx_mult(z, z);

  z1 = cmplx_sclr(0.5641896, u);
  z1 = cmplx_mult(z, cmplx_addr(z1, 1.410474));

  z2 = cmplx_mult(u, cmplx_addr(u, 3.0));
  z2 = cmplx_addr(z2, 0.75);

  return cmplx_div(z1, z2);
}
/* ------- end ---------------------------- Humlicek2.c ------------- */

#define N_HUMLICEK_3  5

/* ------- begin -------------------------- Humlicek3.c ------------- */

complex Humlicek3(complex z)
{
  register int n;

  static double a[N_HUMLICEK_3] =
                     {0.5642236, 3.778987, 11.96482, 20.20933, 16.4955};
  static double b[N_HUMLICEK_3] =
                     {6.699398, 21.69274,  39.27121, 38.82363, 16.4955};

  complex z1, z2;

  /* --- Approximation in region III --                -------------- */

  z1 = cmplx(a[0], 0.0);
  z2 = cmplx_addr(z, b[0]);

  for (n = 1;  n < N_HUMLICEK_3;  n++) {
    z1 = cmplx_addr(cmplx_mult(z1, z), a[n]);
    z2 = cmplx_addr(cmplx_mult(z2, z), b[n]);
  }

  return cmplx_div(z1, z2);
}
/* ------- end ---------------------------- Humlicek3.c ------------- */

#define N_HUMLICEK_4  7

/* ------- begin -------------------------- Humlicek4.c ------------- */

complex Humlicek4(complex z)
{
  register int n;

  static double a[N_HUMLICEK_4] =
     {0.56419, 1.320522, 35.7668, 219.031, 1540.787, 3321.99, 36183.31};
  static double b[N_HUMLICEK_4] =
     {1.841439, 61.57037, 364.2191, 2186.181, 9022.228, 24322.84, 32066.6};

  complex z1, z2, u, expu;

  u = cmplx_mult(z, cmplx_sclr(-1.0, z));

  /* --- Approximation in region IV --                 -------------- */

  z1 = cmplx(a[0], 0.0);
  z2 = cmplx_addr(u, b[0]);

  for (n = 1;  n < N_HUMLICEK_4;  n++) {
    z1 = cmplx_addr(cmplx_mult(u, z1), a[n]);
    z2 = cmplx_addr(cmplx_mult(u, z2), b[n]);
  }
  expu = cmplx_exp(cmplx_sclr(-1.0, u));
  
  return cmplx_subt(expu, cmplx_div(cmplx_mult(z, z1), z2));
}
/* ------- end ---------------------------- Humlicek4.c ------------- */
