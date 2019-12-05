/* ------- file: -------------------------- voigt.c -----------------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Fri Dec 14 08:53:28 2001 --

       --------------------------                      ----------RH-- */

/* --- Voigt profile generators:

       -- Armstrong 1967, JQSRT 7, pp. 61-88
          (slow for damping parameters larger than 1.5, accurate to
          6 significant figures).

       -- Hui, Armstrong & Wray 1978, JQSRT 19, pp. 509-516
          (same speed in whole parameter space, faster than Armstrong for
          large a, otherwise comparable, only 1% accurate for larger v).

       -- George Rybicki's accurate (Ref)
          (accurate to at least 8 significant figures, but a factor 2
          slower than Armstrong and Hui et al.).

       -- David Hummer's (Ref: algorithm 363 from the collected algorithms
          from CACM, and SIAM J. Num. Analysis, 7, 185 (1970))
          (Speed is comparable to Rybicki's).

       -- Humlicek 1982, JQSRT 27, p. 437
          Relative accuracy 1.0E-04. Also calculates Faraday-Voigt
          function needed in Stokes radiative transfer.

       -- Lookup table. 

    Note: If a FORTRAN 90 compiler is available FORTRAN's builtin
          complex arithmatic can be used by defining HAVE_F90 and
          linking with humlicek_.f90, hui_.f90, and libf90 at the final
          linking step that makes the executable.

       --                                              -------------- */

 
#include <stdlib.h>
#include <math.h>

#include "rh.h"
#include "constant.h"
#include "complex.h"
#include "error.h"

#define  TINY 1.0E-08


/* --- Function prototypes --                          -------------- */

double VoigtArmstrong(double a, double v);
double VoigtRybicki(double a, double v);
double VoigtHui(double a, double v, double *F);
double VoigtHumlicek(double a, double v, double *F);

double VoigtK1(double a, double v);
double VoigtK2(double a, double v);
double VoigtK3(double a, double v);

#if defined(HAVE_F90)
void humlicek_(double *a, double *v, complex *W);
void hui_(complex *z, complex *W);
#else
complex Humlicek1(complex z);
complex Humlicek2(complex z);
complex Humlicek3(complex z);
complex Humlicek4(complex z);
#endif

double VoigtLookup(double a, double v);


/* --- Global variables --                             -------------- */

extern char messageStr[];


/* ------- begin -------------------------- Voigt.c ----------------- */

double Voigt(double a, double v, double *F, enum VoigtAlgorithm algorithm)
{
  double voigt;

  switch (algorithm) {
  case ARMSTRONG:
    voigt = VoigtArmstrong(a, v);
    break;
  case RYBICKI:
    voigt = VoigtRybicki(a, v);
    break;
  case HUI_ETAL:
    voigt = VoigtHui(a, v, F);
    break;
  case HUMLICEK:
    voigt = VoigtHumlicek(a, v, F);
    break;
  case LOOKUP:
    voigt = VoigtLookup(a, v);
    break;
  default:
    sprintf(messageStr, "Unregognized Voigt algorithm: %d", algorithm);
    Error(ERROR_LEVEL_2, "Voigt", messageStr);
  }
  return voigt;
}
/* ------- end ---------------------------- Voigt.c ----------------- */

/* ------- begin -------------------------- VoigtArmstrong.c -------- */

#define  NT  10

static double t[NT] = 
               {0.2453407083, 0.7374737285, 1.2340762153, 1.7385377121,
	        2.2549740020, 2.7888060584, 3.3478545673, 3.9447640401,
	        4.6036824495, 5.3874808900};

static double w[NT] =
               {4.6224366960e-01, 2.8667550536e-01, 1.0901720602e-01,
		2.4810520887e-02, 3.2437733422e-03, 2.2833863601e-04,
		7.8025564785e-06, 1.0860693707e-07, 4.3993409922e-10,
		2.2293936455e-13};

double VoigtArmstrong(double a, double v)
{
  if (v < 0.0) v = -v;

  if ((a < 1.0 && v < 4.0)  ||  (a < 1.8/(v + 1.0)))
    return VoigtK1(a, v);
  else if (a < 2.5  &&  v < 4.0)
    return VoigtK2(a, v);
  else
    return VoigtK3(a, v);
}
/* ------- end ---------------------------- VoigtArmstrong.c -------- */

/* ------- begin -------------------------- VoigtK1.c --------------- */

#define TWOSQRTPI  1.12837917
#define EXPMAX     70.0
#define NC         34

double VoigtK1(double a, double v)
{
  static double c[NC] =
                 { 0.1999999999972224, -0.1840000000029998,
		   0.1558399999965025, -0.1216640000043988,
		   0.0877081599940391, -0.0585141248086907,
		   0.0362157301623914, -0.0208497654398036,
		   0.0111960116346270, -0.56231896167109e-02,
		   0.26487634172265e-02, -0.11732670757704e-02,
		   0.4899519978088e-03, -0.1933630801528e-03,
		   0.722877446788e-04, -0.256555124979e-04,
		   0.86620736841e-05, -0.27876379719e-05,
		   0.8566873627e-06, -0.2518433784e-06, 0.709360221e-07,
		  -0.191732257e-07, 0.49801256e-08, -0.12447734e-08,
		   0.2997777e-09, -0.696450e-10, 0.156262e-10, -0.33897e-11,
		   0.7116e-12, -0.1447e-12, 0.285e-13, -0.55e-14,
		   0.10e-14, -0.2e-15 };

  register int n;

  double a2 = a*a, v2 = v*v, u1, dn01, dn02, dn, v2i, funct, an, q, g,
         coef, bn01, bn02, bn, v1;

  if ((v2 - a2) > EXPMAX)
    u1 = 0.0;
  else
    u1 = exp(a2 - v2) * cos(2.0*v*a);

  if (v > 5.0) {
    v2i = 1.0 / v2;
    dn01 = -v2i * (0.5 + v2i*(0.75 + v2i*(1.875 + v2i*(6.5625 +
           v2i*(29.53125 + v2i*(1162.4218 + v2i*1055.7421))))));
    dn02 = (1.0 - dn01) / (2.0 * v);
  } else {
    bn01 = bn02 = 0.0;
    v1   = v / 5.0;
    coef = 4.0 * v1*v1 - 2.0;
    for (n = NC-1;  n >= 0;  n--) {
      bn = coef*bn01 - bn02 + c[n];
      bn02 = bn01;
      bn01 = bn;
    }
    dn02 = (double) (v1*(bn - bn02));
    dn01 = 1.0 - 2.0*v*dn02;
  }

  funct = a*dn01;
  if (a > TINY) {
    q  = 1.0;
    an = a;
    for (n = 2;  n <= 50;  n++) {
      dn   = (v*dn01 + dn02) * (-2.0/n);
      dn02 = dn01;
      dn01 = dn;
      if (n % 2) {
	q   = -q;
	an *= a2;
	g   = dn * an;
	funct += q*g;
	if (fabs(g/funct) <= TINY) return (u1 - TWOSQRTPI*funct);
      }
    }
  }
  return (u1 - TWOSQRTPI*funct);
}
/* ------- end ---------------------------- VoigtK1.c --------------- */

/* ------- begin -------------------------- VoigtK2.c --------------- */

double VoigtK2(double a, double v)
{
  register int n;

  double g = 0.0, r, s, a2 = a*a;

  for (n = 0;  n < NT;  n++) {
    r  = t[n] - v;
    s  = t[n] + v;
    g += (4.0*t[n]*t[n] - 2.0) * (r*atan(r/a) + s*atan(s/a) - 
                      0.5*a*(log(a2 + r*r) + log(a2 + s*s))) * w[n];
  }
  return  g/PI;
}
/* ------- end ---------------------------- VoigtK2.c --------------- */

/* ------- begin -------------------------- VoigtK3.c --------------- */

double VoigtK3(double a, double v)
{
  register int n;

  double g = 0.0, a2 = a*a;

  for (n = 0;  n < NT;  n++) {
    g += (1.0/(SQ(v - t[n]) + a2) + 1.0/(SQ(v + t[n]) + a2)) * w[n];
  }
  return  (a*g)/PI;
}
/* ------- end ---------------------------- VoigtK3.c --------------- */

/* ------- begin -------------------------- VoigtRybicki.c ---------- */

#define  NGR       31

#define  THREEPI    9.42477796076938
#define  INVSQRTPI  0.564189583547756
#define  C0         0.0897935610625833
#define  C1        29.608813203268


/* --- Computes the Voigt function for any value of x and
       any positive value of a. Adapted from original by G. Rybicki - */

double VoigtRybicki(double a, double v)
{
  register int m, n;
  static   int initialize = TRUE;
  static   double c[NGR];

  double a1, a2, b1, b2, e, s, t, zi, zr, voigt;

  if (initialize) {
    for (m = -15, n = 0;  n < NGR;  m++, n++)
      c[n] = C0 * exp(-(m * m)/9.0);
    initialize = FALSE;
  }
  /* --- Doppler profile (a = 0.0) --                   ------------- */

  if (a == 0.0)
    return  (double) exp(-v*v);

  /* --- General case --                                ------------- */

  a1 = 3.0 * a;
  a2 = a * a;
  e  = exp(-THREEPI*a);
  if (a < 0.1) {
    zr = 0.5 * (e + 1.0/e) * cos(THREEPI * v);
    zi = 0.5 * (e - 1.0/e) * sin(THREEPI * v);
    voigt = INVSQRTPI * exp(a2 - v*v) * cos(2.0 * a*v);
  } else {
    zr = e * cos(THREEPI * v);
    zi = e * sin(THREEPI * v);
    voigt = 0.0;
  }
  b1 =  (1.0 - zr) * a * 1.5;
  b2 = -zi;
  s  = -8.0 - 1.5*v;
  t  = s*s + 2.25*a2;
  for (n = 0;  n < NGR;  n++) {
    t  =  t + s + 0.25;
    s  =  s + 0.5;
    b1 =  a1 - b1;
    b2 = -b2;
    if (t > 2.5e-12)
      voigt += c[n] * (b1 + b2*s) / t;
    else
      voigt -= c[n] * a * C1;
  }
  return (double) voigt * SQRTPI;
}
/* ------- end ---------------------------- VoigtRybicki.c ---------- */

/* ------- begin -------------------------- VoigtHui.c -------------- */

#define NHUI 6

double VoigtHui(double a, double v, double *F)
{
#if !defined(HAVE_F90)
  register int n;

  static double ah[NHUI+1] = 
        {122.607931777104326, 214.382388694706425, 181.928533092181549,
          93.155580458138441,  30.180142196210589,   5.912626209773153,
           0.564189583562615};

  static double bh[NHUI+1] =
        {122.607931773875350, 352.730625110963558, 457.334478783897737,
         348.703917719495792, 170.354001821091472,  53.992906912940207,
          10.479857114260399};

  complex W1 = {0.0, 0.0}, W2;

  /* --- Voigt function generator using rational approximation of
         the complex error function W(a - iv).
         This routine is faster for large a (> 1.5)

         Voigt:             H(a, v) = Re[W(v + ia)]
         Faraday - Voigt: 2*F(a, v) = Im[W(v + ia)]

         --                                            -------------- */
#endif
  complex z, W;

  z  = cmplx(a, -v);

#if defined(HAVE_F90)

  /* --- If a FORTRAN 90 compiler is available then FORTRAN's intrinsic 
         complex arithmatic can be used. This will give a factor of
         two to three improvement in speed (with the SUN compilers).

    See: hui_.f90
         --                                            -------------- */

  hui_(&z, &W);
#else

  /* --- C version is used otherwise --                -------------- */

  W2 = z;
  for (n = NHUI;  n >= 0;  n--) {
    W1 = cmplx_mult(cmplx_addr(W1, ah[n]), z);
    W2 = cmplx_mult(cmplx_addr(W2, bh[n]), z);
  }
  W = cmplx_div(W1, W2);
#endif

  if (F != NULL) *F = W.i;
  return W.r;
}
/* ------- end ---------------------------- VoigtHui.c -------------- */

/* ------- begin -------------------------- VoigtHumlicek.c --------- */

  /* --- Voigt function generator using rational approximation of
         the complex error function W(v + ia).

         Output:

         Voigt:             H(a, v) = Re[W(v + ia)]
         Faraday - Voigt: 2*F(a, v) = Im[W(v + ia)]

         --                                            -------------- */

double VoigtHumlicek(double a, double v, double *F)
{
#if defined(HAVE_F90)

  complex W;

  /* --- If a FORTRAN 90 compiler is available then FORTRAN's intrinsic 
         complex arithmatic can be used. This will give a factor of
         two to three improvement in speed (with the SUN compilers).

    See: humlicek_.f90
         --                                            -------------- */

  humlicek_(&a, &v, &W);

#else

  complex W, z = cmplx(a, -v);
  double s = fabs(v) + a;

  /* --- C versions are called otherwise

    See: humlicek.c
         --                                            -------------- */

  if (s >= 15.0)
    W = Humlicek1(z);
  else if (s >= 5.5)
      W = Humlicek2(z);
  else if (a >= 0.195*fabs(v) - 0.176)
    W = Humlicek3(z);
  else
    W = Humlicek4(z);

#endif

  if (F != NULL) *F = W.i;
  return W.r;
}
/* ------- end ---------------------------- VoigtHumlicek.c --------- */

/* ------- begin -------------------------- VoigtLookup.c ----------- */

#define TABLE_ALGORITHM  RYBICKI

#define A_MIN  1.0E-05
#define A_MAX  1.0E-02
#define N_A    50

#define V_MIN_LIN    0.0
#define V_MAX_LIN    4.0
#define N_V_LIN      400    
#define V_MIN_LOG V_MAX_LIN
#define V_MAX_LOG 1000.0
#define N_V_LOG       30
     

/* --- Compute lookup table of Voigt function when called
       for the first time. With subsequent calls the lookup table is
       interpolated (bi-linearly for the moment) --    -------------- */

double VoigtLookup(double a, double v)
{
  const char routineName[] = "VoigtLookup";
  register int n, m;

  static bool_t initialize = TRUE, hunt;
  static double *a_table, *v_table_lin, *v_table_log,
               **table_lin, **table_log, voigt, log_a;

  double da, dv;

  v = fabs(v);

  if (a < A_MIN || a > A_MAX || v > V_MAX_LOG) {
    sprintf(messageStr,
	    "Arguments (%E, %E) outside domain table", a, v);
    Error(WARNING, routineName, messageStr);

    return Voigt(a, v, NULL, TABLE_ALGORITHM);
  }

  /* --- Populate the lookup tables on initialization -- ------------ */

  if (initialize) {
    a_table     = (double *) malloc(N_A * sizeof(double));
    v_table_lin = (double *) malloc(N_V_LIN * sizeof(double));
    v_table_log = (double *) malloc(N_V_LOG * sizeof(double));

    table_lin = matrix_double(N_A, N_V_LIN);
    table_log = matrix_double(N_A, N_V_LOG);

    a_table[0]     = log(A_MIN);
    a_table[N_A-1] = log(A_MAX);
    da = (a_table[N_A-1] - a_table[0]) / (N_A - 1);
    for (n = 1;  n < N_A-1;  n++)
      a_table[n] = a_table[n-1] + da;

    v_table_lin[0]         = V_MIN_LIN;
    v_table_lin[N_V_LIN-1] = V_MAX_LIN;
    dv = (v_table_lin[N_V_LIN-1] - v_table_lin[0]) / (N_V_LIN - 1);
    for (m = 1;  m < N_V_LIN-1;  m++)
      v_table_lin[m] = v_table_lin[m-1] + dv;

    v_table_log[0]         = log(V_MIN_LOG);
    v_table_log[N_V_LOG-1] = log(V_MAX_LOG);
    dv = (v_table_log[N_V_LOG-1] - v_table_log[0]) / (N_V_LOG - 1);
    for (m = 1;  m < N_V_LOG-1;  m++)
      v_table_log[m] = v_table_log[m-1] + dv;

    for (n = 0;  n < N_A;  n++) {
      for (m = 0;  m < N_V_LIN;  m++) {
        table_lin[n][m] = Voigt(exp(a_table[n]), v_table_lin[m],
				NULL, TABLE_ALGORITHM);
      }
      for (m = 0;  m < N_V_LOG;  m++) {
        table_log[n][m] = Voigt(exp(a_table[n]), exp(v_table_log[m]),
				NULL, TABLE_ALGORITHM);
      }
    }
    sprintf(messageStr, "Created Voigt lookup tables\n");
    Error(MESSAGE, routineName, messageStr);

    initialize = FALSE;
  }
  /* --- There are two domains. Interpolation is always logarithmic
         in damping parameter a. For large values of Doppler frequency
         v it is also logarithmic in v, otherwise it is linear in v - */

  log_a = log(a);

  if (v < V_MAX_LIN) {
    voigt = BiLinear(N_A, a_table, log_a, N_V_LIN, v_table_lin, v,
		     table_lin, hunt=TRUE);
  } else {
    voigt = BiLinear(N_A, a_table, log_a, N_V_LOG, v_table_log, log(v),
		     table_log, hunt=FALSE);
  }
  return voigt;
}
/* ------- begin -------------------------- VoigtLookup.c ----------- */
