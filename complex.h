/* ------- file: -------------------------- complex.h ---------------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Mon Nov 22 15:35:13 1999 --

       --------------------------                      ----------RH-- */

#ifndef __COMPLEX_H__
#define __COMPLEX_H__

/* --- Defines structure for complex numbers.
       Let z = a + b*i, then the real part is represented as z.r, 
       the imaginary part as z.i. Initialize as cmplx z = {a, b}. --  */

typedef struct {double r; double i;} complex; 

#define FFT_FORWARD -1
#define FFT_REVERSE  1
#define FFT_REAL     0
#define FFT_COMPLEX  1


/* --- Associated function prototypes --               -------------- */

complex **matrix_complex(int Nrow, int Ncol);


/* --- Manipulation of complex numbers --              -------------- */

complex cmplx(double a, double b);
complex cmplx_conj(complex a);
complex cmplx_mult(complex a, complex b);
complex cmplx_sclr(double a, complex b);
complex cmplx_div(complex a, complex b);
complex cmplx_exp(complex a);
complex cmplx_add(complex a, complex b);
complex cmplx_addr(complex a, double b);
complex cmplx_subt(complex a, complex b);


/* --- Fourier transforms --                           -------------- */

void  four1(int NN, double *data, int sign);
void  fourn(int Ndim, int *NN, double *data, int sign);

void  fourier_(int *Ndim, int *NN, double *data, int *sign, 
               int *format, double *work);

#endif /* !__COMPLEX_H__ */

/* ------- end ---------------------------- complex.h --------------- */
