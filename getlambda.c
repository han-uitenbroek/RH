/* ------- file: -------------------------- getlambda.c -------------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Wed Apr 22 08:59:40 2009 --

       --------------------------                      ----------RH-- */

/* --- Construct a wavelength grid that is approximately equidistant
       in the core (q <= qcore) and equidistant in log(q) in the wings
       (qcore < q <= qwing).

       Look for a function of the form: q[n] = a*(n + (exp(b*n)-1));
       n=0, N-1 satisfying the following conditions:

    -- q[0]   = 0          (this is true for all a and b)
    -- q[(N-1)/2] = qcore
    -- q[N-1] = qwing

    --                                                 -------------- */
 
#include <math.h>
#include <stdlib.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "constant.h"
#include "error.h"
#include "inputs.h"

/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern InputData input;
extern char   messageStr[];


/* ------- begin -------------------------- getLambda.c ------------- */

void getLambda(AtomicLine *line)
{
  const char routineName[] = "getLambda";
  register int la, n;

  int     Nlambda, NB = 0, Nmid;
  double  beta, a, b, y, q_to_lambda, *q, qB_char, qB_shift, dlambda,
          g_Lande_eff, lambda0;

  /* --- First some sanity checks on qcore and qwing -- ------------- */

  if ((line->qcore <= 0.0)  ||  (line->qwing <= 0.0)){
    sprintf(messageStr, "Either qcore or qwing is negative or zero"
                        " transition %d->%d", line->j, line->i);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }
  /* --- Create full (instead of half) line profiles in case of:

         1) a set keyword ASYMM in this_atom.atom for this line
	 2) a moving atmosphere
	 3) a magnetic field is present and the g_Lande is
             unequal 0.0 for this line
         4) the line has more than one components

         --                                            -------------- */

  if ((atmos.Stokes && line->polarizable) ||
      atmos.moving ||
      line->Ncomponent > 1) {
    line->symmetric = FALSE;
  }
  /* --- Make sure the half-line profile always has an odd number of
         grid points --                                -------------- */

  if (line->symmetric) {
    Nlambda = (line->Nlambda % 2) ? line->Nlambda : line->Nlambda + 1;
  } else {
    Nlambda =
      (line->Nlambda % 2) ? (line->Nlambda / 2) : (line->Nlambda + 1)/2;
  }

  if (line->qwing <= 2.0*line->qcore) {

    /* --- In this case use a linear scale out to qwing -- ---------- */

    sprintf(messageStr, "Ratio of qwing / (2*qcore) is less than one.\n"
                        " Using linear spacing for transition: %d->%d",
	    line->j, line->i);
    Error(WARNING, routineName, messageStr);
    beta = 1.0;
  } else
    beta = line->qwing / (2.0*line->qcore);

  y = beta + sqrt(SQ(beta) + (beta - 1.0)*Nlambda + 2.0 - 3.0*beta);
  b = 2.0*log(y) / (Nlambda - 1);
  a = line->qwing / (Nlambda - 2.0 + SQ(y));

  q = (double *) malloc(Nlambda * sizeof(double));
  for (la = 0;  la < Nlambda;  la++)
    q[la] = a * (la + (exp(b * la) - 1.0));

  if (line->polarizable) {
    g_Lande_eff = effectiveLande(line);

    /* --- When magnetic field is present account for denser
           wavelength spacing in the unshifted \pi component, and the
           red and blue-shifted circularly polarized components.
           First, get characteristic Zeeman splitting -- ------------ */  
    
    qB_char = g_Lande_eff * (Q_ELECTRON / (4.0*PI * M_ELECTRON)) *
      (line->lambda0 * NM_TO_M) * atmos.B_char / atmos.vmicro_char;

    if (qB_char/2.0 >= line->qcore) {
      sprintf(messageStr,
      "Characteristic Zeeman splitting qB_char (= %f) >= 2.0*qcore for "
      "transition %d->%d", qB_char, line->j, line->i);
      Error(WARNING, routineName, messageStr);
    }
    Locate(Nlambda, q, qB_char/2.0, &NB);
    qB_shift = 2 * q[NB];

    q = realloc(q, (Nlambda + 2*NB) * sizeof(double));
    for (la = NB+1;  la <= 2*NB;  la++)
      q[la] = qB_shift - a * (2*NB-la + (exp(b * (2*NB-la)) - 1.0));

    for (la = 2*NB+1;  la < Nlambda + 2*NB;  la++)
      q[la] = qB_shift + a * (la - 2*NB + (exp(b * (la - 2*NB)) - 1.0));

    Nlambda += 2*NB;
  }
  /* --- Store absolute wavelength rather than differences -- ----- */

  q_to_lambda = line->lambda0 * (atmos.vmicro_char / CLIGHT);

  if (line->symmetric) {
    line->Nlambda = Nlambda;
    line->lambda = (double *) malloc(line->Nlambda * sizeof(double));
    for (la = 0;  la < line->Nlambda;  la++)
      line->lambda[la] = line->lambda0 + q_to_lambda * q[la];
  } else {
    line->Nlambda = (2*Nlambda - 1) * line->Ncomponent;
    line->lambda = (double *) malloc(line->Nlambda * sizeof(double));

    for (n = 0;  n < line->Ncomponent;  n++) {
      Nmid = n*(2*Nlambda - 1) + Nlambda - 1;
      lambda0 = line->lambda0 + line->c_shift[n];

      line->lambda[Nmid] = lambda0;
      for (la = 1;  la < Nlambda;  la++) {
	dlambda = q_to_lambda * q[la];
	line->lambda[Nmid - la] = lambda0 - dlambda;
	line->lambda[Nmid + la] = lambda0 + dlambda;
      }
    }
    if (line->Ncomponent > 1)
      qsort(line->lambda, line->Nlambda, sizeof(double), qsascend);
  }
  free(q);
}
/* ------- end ---------------------------- getLambda.c ------------- */

/* ------- begin -------------------------- getLambdaCont.c --------- */

void getLambdaCont(AtomicContinuum *continuum, double lambdamin)
{
  register int la;

  int    Nlambda;
  double dlamb;

  Nlambda = continuum->Nlambda;
  dlamb   = (continuum->lambda0 - lambdamin) / (Nlambda - 1);

  continuum->lambda[0] = lambdamin;
  for (la = 1;  la < Nlambda;  la++)
    continuum->lambda[la] = continuum->lambda[la-1] + dlamb;
}
/* ------- end ---------------------------- getLambdaCont.c --------- */

/* ------- begin -------------------------- getwlambda_line.c ------- */

double getwlambda_line(AtomicLine *line, int la)
{
  double wlambda = 0.0, Dopplerwidth;

  /* --- Return wavelength interval. In case of atomic bound-bound
         or molecular vibration-rotation transition return interval
         in Doppler units (without factoring in the thermal velocity  */

  Dopplerwidth = CLIGHT / line->lambda0;
  if (line->symmetric) Dopplerwidth *= 2.0;

  if (la == 0)
    wlambda = 0.5 * (line->lambda[la+1] - line->lambda[la]);
  else if (la == line->Nlambda-1)
    wlambda = 0.5 * (line->lambda[la] - line->lambda[la-1]);
  else
    wlambda = 0.5 * (line->lambda[la+1] - line->lambda[la-1]);

  return wlambda * Dopplerwidth;
}
/* ------- end ---------------------------- getwlambda_line.c ------- */

/* ------- begin -------------------------- getwlambda_cont.c ------- */

double getwlambda_cont(AtomicContinuum *continuum, int la)
{
  double wlambda = 0.0;

  if (la == 0)
    wlambda = 0.5 * (continuum->lambda[la+1] - continuum->lambda[la]);
  else if (la == continuum->Nlambda-1)
    wlambda = 0.5 * (continuum->lambda[la] - continuum->lambda[la-1]);
  else
    wlambda = 0.5 * (continuum->lambda[la+1] - continuum->lambda[la-1]);

  return wlambda;
}
/* ------- end ---------------------------- getwlambda_cont.c ------- */

/* ------- begin -------------------------- getwlambda_mrt.c -------- */

double getwlambda_mrt(MolecularLine *mrt, int la)
{
  double wlambda = 0.0, Dopplerwidth;

  Dopplerwidth = CLIGHT / mrt->lambda0;
  if (mrt->symmetric) Dopplerwidth *= 2.0;

  if (la == 0)
    wlambda = 0.5 * (mrt->lambda[la+1] - mrt->lambda[la]);
  else if (la == mrt->Nlambda-1)
    wlambda = 0.5 * (mrt->lambda[la] - mrt->lambda[la-1]);
  else
    wlambda = 0.5 * (mrt->lambda[la+1] - mrt->lambda[la-1]);

  return wlambda * Dopplerwidth;
}
/* ------- end ---------------------------- getwlambda_mrt.c -------- */
