/* ------- file: -------------------------- cubeconvol.c ------------

       Version:       rh2.0
       Author:        Han Uitenbroek  (huitenbroek@nso.edu)
       Last modified: Mon Sep 12 22:48:11 2005 --

       --------------------------                      ----------RH-- */

/* --- Routines for interpolation by cubic convolution.

  See: R.G. Keys, 1981, in IEEE Trans. Acoustics, Speech,
        and Signal Processing, Vol. 29, pp. 1153-1160.
       --                                              -------------- */
 
#include <math.h>
#include <stdlib.h>

#include "rh.h"

#define  NCC 4


/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */


/* ------- begin -------------------------- cubeconvol.c ------------ */

double cubeconvol(int Nx, int Ny, double *f, double x, double y)
{
  register int n, m;

  int    i, j;
  double *fc, ip, jp, ux[NCC], uy[NCC], c[NCC][NCC], g;

  if (x <= 0.0) {
    i = 0;
    cc_kernel(0.0, ux);
  } else if (x >= (double) (Nx - 1)) {
    i = Nx - 2;
    cc_kernel(1.0, ux);
  } else {
    cc_kernel(modf(x, &ip), ux);
    i = (int) ip;
  }
  if (y <= 0.0) {
    j = 0;
    cc_kernel(0.0, uy);
  } else if (y >= (double) (Ny - 1)) {
    j = Ny - 2;
    cc_kernel(1.0, uy);
  } else {
    cc_kernel(modf(y, &jp), uy);
    j = (int) jp;
  }
 
  g  = 0.0;
  fc = f + (Nx*(j-1) + i-1);

  if (j == 0) {

    /* --- Interpolation point at lower boundary --    -------------- */

    fc += Nx;
    if (i == 0) {
      for (m = 1;  m < NCC;  m++, fc += Nx) {
	for (n = 1;  n < NCC;  n++) c[m][n] = fc[n];
	c[m][0] = 3.0*(c[m][1] - c[m][2]) + c[m][3];
      }
    } else if (i == Nx-2) {
      for (m = 1;  m < NCC;  m++, fc += Nx) {
	for (n = 0;  n < NCC-1;  n++) c[m][n] = fc[n];
	c[m][3] = 3.0*(c[m][2] - c[m][1]) + c[m][0];
      }
    } else {
      for (m = 1;  m < NCC;  m++, fc += Nx)
	for (n = 0;  n < NCC;  n++) c[m][n] = fc[n];
    }
    for (n = 0;  n < NCC;  n++)
      c[0][n] = 3.0*(c[1][n] - c[2][n]) + c[3][n];
  } else if (j == Ny-2) {

    /* --- Upper boundary --                          --------------- */

    if (i == 0) {
      for (m = 0;  m < NCC-1;  m++, fc += Nx) {
	for (n = 1;  n < NCC;  n++) c[m][n] = fc[n];
	c[m][0] = 3.0*(c[m][1] - c[m][2]) + c[m][3];
      }
    } else if (i == Nx-2) {
      for (m = 0;  m < NCC-1;  m++, fc += Nx) {
	for (n = 0;  n < NCC-1;  n++) c[m][n] = fc[n];
	c[m][3] = 3.0*(c[m][2] - c[m][1]) + c[m][0];
      }
    } else {
      for (m = 0;  m < NCC-1;  m++, fc += Nx)
	for (n = 0;  n < NCC;  n++) c[m][n] = fc[n];
    }
    for (n = 0;  n < NCC;  n++)
      c[3][n] = 3.0*(c[2][n] - c[1][n]) + c[0][n];
  } else {

    /* --- General case in vertical direction --      --------------- */

    if (i == 0) {
      for (m = 0;  m < NCC;  m++, fc += Nx) {
	for (n = 1;  n < NCC;  n++) c[m][n] = fc[n];
	c[m][0] = 3.0*(c[m][1] - c[m][2]) + c[m][3];
      }
    } else if (i == Nx-2) {
      for (m = 0;  m < NCC;  m++, fc += Nx) {
	for (n = 0;  n < NCC-1;  n++) c[m][n] = fc[n];
	c[m][3] = 3.0*(c[m][2] - c[m][1]) + c[m][0];
      }
    } else {

      /* --- When interpolation point lies in interior then convolute
             without copying into c, and return immediately -- ------ */

      for (m = 0;  m < NCC;  m++, fc += Nx) {
	for (n = 0;  n < NCC;  n++) {
	  g += fc[n] * ux[n] * uy[m];
	}
      }
      return g;
    }
  }
  /* --- Do the convolution for the boundary cases --  -------------- */

  for (m = 0;  m < NCC;  m++) {
    for (n = 0;  n < NCC;  n++) {
      g += c[m][n] * ux[n] * uy[m];
    }
  }
  return g;
}
/* ------- end ---------------------------- cubeconvol.c ------------ */

/* ------- begin -------------------------- cc_kernel.c ------------- */

void cc_kernel(double s, double *kernel)
{
  /* --- Interpolation kernel for cubic convolution --  ------------- */

  kernel[0] = -s * (1.0 - s*(2.0 - s)) / 2.0;
  kernel[1] = 1.0 + s*s * (3.0*s - 5.0) / 2.0;
  kernel[2] = s * (1.0 + s*(4.0 - 3.0*s)) / 2.0;
  kernel[3] = s*s * (s - 1.0) / 2.0;
}
/* ------- end ---------------------------- cc_kernel.c ------------- */
