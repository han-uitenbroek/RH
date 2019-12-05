/* ------- file: -------------------------- giigen.c ----------------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Thu Jul  1 05:11:13 2010 --

       --------------------------                      ----------RH-- */

/* --- Evaluates angle-averaged general redistribution function PII for
       ordinary and cross redistribution of lines with sharp lower and
       broadened upperlevels.

  See: D. G. Hummer 1962, MNRAS, 125, 21-37.
       I. Huben\'y, 1982, JQSRT 27, 593.

       Parameters:

         adamp     --  Damping parameter line.
         waveRatio --  Ratio of wavelengths of lines in which absorption and
                       emission occurs respectively (set to 1.0 for ordinary
                       redistribution).
         q_abs     --  Frequency absorbed photon (Doppler units)
         q_emit    --  Frequency emitted photon (Doppler units) -- --- */

#include <math.h>
#include <stdio.h>

#include "rh.h"
#include "atom.h"
#include "constant.h"

/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */


/* ------- begin -------------------------- PII --------------------- */

#define NGAUSS   8

double PII(double adamp, double waveRatio, double q_emit, double q_abs)
{
  register int  n;

  static bool_t initialize = TRUE;
  double theta, pii, a1, a2, b1, b2, c1;

  /* --- Use 8-point Gaussian quadrature --           --------------- */

  static double xg[NGAUSS/2] =
    {0.183434642495, 0.525532409916, 0.796666477413, 0.960289856497};
  static double wg[NGAUSS/2] =
    {0.362683783378, 0.313706645877, 0.222381034453, 0.101228536290};
  static double sn[NGAUSS], cs[NGAUSS];

  if (initialize) {
    for (n = 0;  n < NGAUSS;  n++) {
      theta = (n % 2)  ? 0.5*PI*(1.0 - xg[n/2]) : 0.5*PI*(1.0 + xg[n/2]);
      sn[n] = sin(theta);
      cs[n] = cos(theta);
    }
    initialize = FALSE;
  }

  pii = 0.0;
  for (n = 0;  n < NGAUSS;  n++) {
    a1 = 1.0 / sqrt(1.0 - waveRatio*(2.0*cs[n] - waveRatio));
    a2 = a1 * waveRatio;
    b1 = (a2 - a1*cs[n]) / sn[n];
    b2 = (a1 - a2*cs[n]) / sn[n];
    c1 = a1*q_abs - a2*q_emit;

    pii += Voigt(adamp/(a2 * sn[n]), b1*q_abs + b2*q_emit, NULL, RYBICKI) *
      exp(-c1*c1) * wg[n/2];
  }
  return (0.25 * pii);
}
/* ------- end ---------------------------- PII --------------------- */

/* ------- begin -------------------------- GII --------------------- */
/*
 * Gouttebroze's fast approximation for
 *  GII(q_abs, q_emit) = PII(q_abs, q_emit) / phi(q_emit)

 * See: P. Gouttebroze, 1986, A&A 160, 195
 *      H. Uitenbroek,  1989, A&A, 216, 310-314 (cross redistribution)
 */

#define GZERO(x)  (1.0 / (fabs(x) + sqrt((x)*(x) + 1.273239545)))

double GII(double adamp, double waveratio, double q_emit, double q_abs)
{
  double gii, pcore, aq_emit, umin, epsilon, giiwing, u1, phicore,
         phiwing; 

  /* --- Symmetrize with respect to emission frequency --   --------- */

  if (q_emit < 0.0) {
    q_emit = -q_emit;
    q_abs  = -q_abs;
  }
  pcore = 0.0;
  gii = 0.0;

  /* --- Core region --                                     --------- */

  if (q_emit < PRD_QWING) {
    if ((q_abs < -PRD_QWING) ||
	(q_abs > q_emit + waveratio*PRD_QSPREAD)) return gii;
    if (fabs(q_abs) <= q_emit)
      gii = GZERO(q_emit);
    else
      gii = exp(SQ(q_emit) - SQ(q_abs)) * GZERO(q_abs);
  
    if (q_emit >= PRD_QCORE) {
      phicore = exp(-SQ(q_emit));
      phiwing = adamp / (SQRTPI * (SQ(adamp) + SQ(q_emit)));
      pcore   = phicore / (phicore + phiwing);
    }
  }
  /* --- Wing region --                                     --------- */

  if (q_emit >= PRD_QCORE) {
    aq_emit = waveratio * q_emit;
    if (q_emit >= PRD_QWING) {
      if (fabs(q_abs - aq_emit) > waveratio*PRD_QSPREAD) return gii;
      pcore = 0.0;
    }
    umin = fabs((q_abs - aq_emit) / (1.0 + waveratio));
    giiwing = (1.0 + waveratio) * (1.0 - 2.0*umin*GZERO(umin)) * 
      exp(-SQ(umin));

    if (waveratio == 1.0) {
      epsilon  = q_abs / aq_emit;
      giiwing *= (2.75 - (2.5 - 0.75*epsilon) * epsilon);
    } else {
      u1 = fabs((q_abs - aq_emit) / (waveratio - 1.0));
      giiwing -= fabs(1.0 - waveratio) * (1.0 - 2.0*u1*GZERO(u1)) *
	exp(-SQ(u1));
    }
    /* --- Linear combination of core- and wing contributions ------- */

    giiwing = giiwing / (2.0 * waveratio * SQRTPI);
    gii = pcore*gii + (1.0 - pcore)*giiwing;
  }
  return gii;
}
/* ------- end ---------------------------- GII --------------------- */
