/* ------- file: -------------------------- complex.c ---------------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Thu Dec 13 11:16:58 2001 --

       --------------------------                      ----------RH-- */

/* --- Complex operations needed in Voigt generator and FFT -- ------ */


#include <stdlib.h>
#include <math.h>

#include "rh.h"
#include "complex.h"
 

/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */


/* ------- begin -------------------------- complex.c --------------- */

complex cmplx(double a, double b)
{
  complex c;

  /* --- Return complex number z = a + bi --             ------------ */

  c.r = a;
  c.i = b;

  return c;
}
/* ------- end ---------------------------- complex.c --------------- */

/* ------- begin -------------------------- cmplx_conj.c ------------ */

complex cmplx_conj(complex a)
{
  complex c;

  /* --- Return complex conjugate --                     ------------ */

  c.r =  a.r;
  c.i = -a.i;

  return c;
}
/* ------- end ---------------------------- cmplx_conj.c ------------ */

/* ------- begin -------------------------- cmplx_mult.c ------------ */

complex cmplx_mult(complex a, complex b)
{
  complex c;

  /* --- Multiply two complex numbers --                 ------------ */

  c.r = a.r*b.r - a.i*b.i;
  c.i = a.r*b.i + a.i*b.r;

  return c;
}
/* ------- end ---------------------------- cmplx_mult.c ------------ */

/* ------- begin -------------------------- cmplx_scl.c ------------- */

complex cmplx_sclr(double a, complex b)
{
  complex c;

  /* --- Multiply a complex number with a scalar --      ------------ */

  c.r = a * b.r;
  c.i = a * b.i;

  return c;
}
/* ------- end ---------------------------- cmplx_scl.c ------------- */

/* ------- begin -------------------------- cmplx_div.c ------------- */

complex cmplx_div(complex a, complex b)
{
  complex c;
  double d;

  /* --- Divide two complex numbers --                   ------------ */

  d = b.r*b.r + b.i*b.i;
  c.r = (a.r*b.r + a.i*b.i) / d;
  c.i = (a.i*b.r - a.r*b.i) / d;

  return c;
}
/* ------- end ---------------------------- cmplx_div.c ------------- */

/* ------- begin -------------------------- cmplx_exp.c ------------- */

complex cmplx_exp(complex a)
{
  /* --- Complex exponential function --                 ------------ */

  return cmplx_sclr(exp(a.r), cmplx(cos(a.i), sin(a.i)));
}
/* ------- end ---------------------------- cmplx_div.c ------------- */

/* ------- begin -------------------------- cmplx_add.c ------------- */

complex cmplx_add(complex a, complex b)
{
  complex c;

  /* --- Add two complex numbers --                      ------------ */

  c.r = a.r + b.r;
  c.i = a.i + b.i;

  return c;
}
/* ------- end ---------------------------- cmplx_add.c ------------- */

/* ------- begin -------------------------- cmplx_addr.c ------------ */

complex cmplx_addr(complex a, double b)
{
  complex c;

  /* --- Add real number to complex --                   ------------ */

  c.r = a.r + b;
  c.i = a.i;

  return c;
}
/* ------- end ---------------------------- cmplx_addr.c ------------ */

/* ------- begin -------------------------- cmplx_subt.c ------------ */

complex cmplx_subt(complex a, complex b)
{
  complex c;

  /* --- Subtract two complex numbers --                 ------------ */

  c.r = a.r - b.r;
  c.i = a.i - b.i;
  return c;
}
/* ------- end ---------------------------- cmplx_subt.c ------------ */

/* ------- begin -------------------------- matrix_complex.c -------- */

complex **matrix_complex(int Nrow, int Ncol)
{
  register int i;

  int      typeSize = sizeof(complex), pointerSize = sizeof(complex *);
  complex *theMatrix, **Matrix;

  theMatrix = (complex *)  calloc(Nrow * Ncol, typeSize);
  Matrix    = (complex **) malloc(Nrow * pointerSize);
  for (i = 0;  i < Nrow;  i++, theMatrix += Ncol)
    Matrix[i] = theMatrix;

  return Matrix;
}
/* ------- end ---------------------------- matrix_complex.c -------- */

