/* ------- file: -------------------------- scatter.c ---------------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Mon May 21 15:18:15 2018 --

       --------------------------                      ----------RH-- */

/* --- Evaluate scattering integral for PRD transition pointed to
       by transition pointer PRDline.

       Adapted for Cross-Redistribution by:
          Eliza Miller-Ricci (Middlebury College), Jun 29 2001 


 Note: Scratch files for the redistribution weights are now written
       in a location determined from PRD_FILE_TEMPLATE. The files are
       no longer automatically deleted. Careful, they can be quite big.
       --                                              -------------- */

 
#include <stdlib.h>
#include <math.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "spectrum.h"
#include "inputs.h"
#include "constant.h"
#include "error.h"
#include "statistics.h"

#define TENSION  8.0
#define PRD_FILE_TEMPLATE "PRD_%s_%d-%d.dat"


/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern Spectrum spectrum;
extern InputData input;
extern char messageStr[];


/* ------- begin -------------------------- PRDScatter.c ------------ */

void PRDScatter(AtomicLine *PRDline, enum Interpolation representation)
{
  const char routineName[] = "scatterIntegral";
  register int  la, k, lap, kr, ip, kxrd;

  char    filename[MAX_LINE_SIZE];
  bool_t  hunt, initialize;
  int     Np, Nread, Nwrite, ij, Nsubordinate;
  double  q_emit, q0, qN, *q_abs = NULL, *qp = NULL, *wq = NULL,
         *qpp = NULL, *gii = NULL, *adamp, cDop, gnorm, *J = NULL,
          Jbar, scatInt, *J_k = NULL, *Pj, gamma, waveratio;
  Atom *atom;
  AtomicLine *line, **XRD, *XRDline;
  AtomicContinuum *continuum;

  /* --- This the static case --                       -------------- */ 

  atom = PRDline->atom;

  if (!PRDline->PRD) {
    sprintf(messageStr, "Line %d -> %d of %2s is not a PRD line",
	    PRDline->j, PRDline->i, atom->ID);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }

  getCPU(3, TIME_START, NULL);

  /* --- Open temporary file for storage of redistribution weights when
         called for the first time  --                 -------------- */

  initialize = FALSE;
  if (PRDline->fp_GII == NULL) {
    sprintf(filename,
	    (atom->ID[1] == ' ') ? "PRD_%.1s_%d-%d.dat" : "PRD_%s_%d-%d.dat",
	    atom->ID, PRDline->j, PRDline->i);

    initialize = TRUE;
    if ((PRDline->fp_GII = fopen(filename, "w+")) == NULL) {
      sprintf(messageStr, "Unable to open temporary file %s", filename);
      Error(ERROR_LEVEL_2, routineName, messageStr);
    }
  }
  if (!initialize) rewind(PRDline->fp_GII);
  
  /* --- Set XRD line array --                         -------------- */

  if (input.XRD)
    Nsubordinate = 1 + PRDline->Nxrd;
  else
    Nsubordinate = 1;

  XRD = (AtomicLine **) malloc(Nsubordinate * sizeof(AtomicLine*));
  XRD[0] = PRDline;
  if (input.XRD) {
    for (kxrd = 0;  kxrd < PRDline->Nxrd;  kxrd++) 
      XRD[kxrd + 1] = PRDline->xrd[kxrd];
  }
  /* --- Initialize the emission profile ratio rho --  -------------- */

  for (la = 0; la < PRDline->Nlambda; la++) {
    for (k = 0; k < atmos.Nspace; k++) {   
      PRDline->rho_prd[la][k] = 1.0;
    }
  }
  Pj    = (double *) malloc (atmos.Nspace * sizeof(double));
  adamp = (double *) malloc (atmos.Nspace * sizeof(double)); 
  cDop  = (NM_TO_M * PRDline->lambda0) / (4.0 * PI);

  for (k = 0; k < atmos.Nspace; k++) { 
    adamp[k] = (PRDline->Grad + PRDline->Qelast[k]) * cDop / atom->vbroad[k];

    /* --- Evaluate the total rate Pj out of the line's upper level - */

    Pj[k] = PRDline->Qelast[k];
    for (ip = 0;  ip < atom->Nlevel;  ip++) {
      ij = ip * atom->Nlevel + PRDline->j;
      Pj[k] += atom->C[ij][k];
    }
    for (kr = 0;  kr < atom->Nline;  kr++) {
      line = &atom->line[kr];
      if (line->j == PRDline->j) Pj[k] += line->Rji[k];
      if (line->i == PRDline->j) Pj[k] += line->Rij[k];
    }
    for (kr = 0;  kr < atom->Ncont;  kr++) {
      continuum = &atom->continuum[kr];
      if (continuum->j == PRDline->j) Pj[k] += continuum->Rji[k];
      if (continuum->i == PRDline->j) Pj[k] += continuum->Rij[k];
    }
  }
  /* --- Loop over subordinate lines --                -------------- */

  for (kxrd = 0;  kxrd < Nsubordinate;  kxrd++) {
    XRDline = XRD[kxrd];
    waveratio = XRDline->lambda0 / PRDline->lambda0;

    J_k   = realloc(J_k, XRDline->Nlambda * sizeof(double));
    q_abs = realloc(q_abs, XRDline->Nlambda * sizeof(double));

    /* --- Loop over all spatial locations --          -------------- */

    for (k = 0;  k < atmos.Nspace;  k++) {
      gamma = atom->n[XRDline->i][k] / atom->n[PRDline->j][k] *
	XRDline->Bij / Pj[k];
      Jbar = XRDline->Rij[k] / XRDline->Bij;

      /* --- Get local mean intensity and wavelength in Doppler units */

      for (la = 0;  la < XRDline->Nlambda;  la++) {
	J_k[la]   = spectrum.J[XRDline->Nblue + la][k];
	q_abs[la] = (XRDline->lambda[la] - XRDline->lambda0) * CLIGHT /
	  (XRDline->lambda0 * atom->vbroad[k]);
      }
      switch (representation) {
      case LINEAR:
	break;
      case SPLINE:
	splineCoef(XRDline->Nlambda, q_abs, J_k);
	break;
      case EXP_SPLINE:
	exp_splineCoef(XRDline->Nlambda, q_abs, J_k, TENSION);
	break;
      }
      /* --- Outer wavelength loop over emission wavelengths -- ----- */

      for (la = 0;  la < PRDline->Nlambda;  la++) {
	q_emit = (PRDline->lambda[la] - PRDline->lambda0) * CLIGHT /
	  (PRDline->lambda0 * atom->vbroad[k]);

	/* --- Establish integration limits over absorption wavelength,
	       using only regions where the redistribution function is
	       non-zero. (See also function GII.) --   -------------- */
 
	if (fabs(q_emit) < PRD_QCORE) {
	  q0 = -PRD_QWING;
	  qN =  PRD_QWING;
	} else {
	  if (fabs(q_emit) < PRD_QWING) {
	    if (q_emit > 0.0) {
	      q0 = -PRD_QWING;
	      qN = waveratio * (q_emit + PRD_QSPREAD);
	    } else {
	      q0 = waveratio * (q_emit - PRD_QSPREAD);
	      qN = PRD_QWING;
	    }
	  } else {
	    q0 = waveratio * (q_emit - PRD_QSPREAD);
	    qN = waveratio * (q_emit + PRD_QSPREAD);
	  }
	}
	Np = (int) ((qN - q0) / PRD_DQ) + 1;
	qp = (double *) realloc(qp,  Np * sizeof(double));
	for (lap = 1, qp[0] = q0;  lap < Np;  lap++)
	  qp[lap] = qp[lap - 1] + PRD_DQ;

	/* --- Fold interpolation of J for symmetric lines. -- ------ */

	if (XRDline->symmetric) {
	  qpp = (double *) realloc(qpp,  Np * sizeof(double));
	  for (lap = 0;  lap < Np;  lap++) qpp[lap] = fabs(qp[lap]);
	}

	/* --- Interpolate mean intensity onto fine grid. Choose linear,
	       spline or exponential spline interpolation -- -------- */

	J = (double *) realloc(J,  Np * sizeof(double));
	switch (representation) {
	case LINEAR:
	  Linear(XRDline->Nlambda, q_abs, J_k, Np,
		 (XRDline->symmetric) ? qpp : qp, J, hunt=TRUE);
	  break;
	case SPLINE:
	  splineEval(Np, (XRDline->symmetric) ? qpp : qp, J, hunt=TRUE);
	  break;
	case EXP_SPLINE:
	  exp_splineEval(Np, (XRDline->symmetric) ? qpp : qp, J, hunt=TRUE);
	  break;
	}
	/* --- Compute the redistribution weights --   -------------- */

	gii = (double *) realloc(gii, Np * sizeof(double));
	if (initialize) {

	  /* --- Integration weights (See: Press et al. Numerical Recipes,
	         p. 107, eq. 4.1.12) --                -------------- */

	  wq = (double *) realloc(wq,  Np * sizeof(double));
	  wq[0] = 5.0/12.0  * PRD_DQ;
	  wq[1] = 13.0/12.0 * PRD_DQ;
	  for (lap = 2;  lap < Np-2;  lap++) wq[lap] = PRD_DQ;
	  wq[Np-1] = 5.0/12.0  * PRD_DQ;
	  wq[Np-2] = 13.0/12.0 * PRD_DQ;

	  for (lap = 0;  lap < Np;  lap++)
	    gii[lap] = GII(adamp[k], waveratio, q_emit, qp[lap]) * wq[lap];

	  if ((Nwrite =
	       fwrite(gii, sizeof(double), Np, PRDline->fp_GII)) != Np) {
	    sprintf(messageStr,
		  "Unable to write proper number of redistribution weights\n"
		  " Wrote %d instead of %d.\n Line %d -> %d, la = %d, k = %d",
		    Nwrite, Np, XRDline->j, XRDline->i, la, k);
	    Error(ERROR_LEVEL_2, routineName, messageStr);
	  }
	} else {
	  if ((Nread =
	       fread(gii, sizeof(double), Np, PRDline->fp_GII)) != Np) {
	    sprintf(messageStr,
		  "Unable to read proper number of redistribution weights\n"
		  " Read %d instead of %d.\n Line %d -> %d, la = %d, k = %d",
		  Nread, Np, XRDline->j, XRDline->i, la, k);
	    Error(ERROR_LEVEL_2, routineName, messageStr);
	  }
	}
	/* --- Inner wavelength loop doing actual wavelength integration
               over absorption wavelengths --          -------------- */

	gnorm   = 0.0;
	scatInt = 0.0;
	for (lap = 0;  lap < Np;  lap++) {
	  gnorm   += gii[lap];
	  scatInt += J[lap] * gii[lap];
	}
	PRDline->rho_prd[la][k] += gamma*(scatInt/gnorm - Jbar);
      }
    }
  }
  /* --- Clean temporary variable space --             -------------- */

  for (kxrd = 0;  kxrd < Nsubordinate;  kxrd++) {
    if (XRD[kxrd]->symmetric) {
      free(qpp);
      break;
    }
  }
  free(J_k);    free(q_abs);  free(gii);
  free(qp);     free(wq);     free(J);
  free(XRD);    free(Pj);     free(adamp);

  sprintf(messageStr, "Scatter Int %5.1f", PRDline->lambda0);
  getCPU(3, TIME_POLL, messageStr);
}
/* ------- end ---------------------------- PRDScatter.c ------------ */

/* ------- begin -------------------------- PRDAngleScatter.c ------- */

void PRDAngleScatter(AtomicLine *PRDline,
		     enum Interpolation representation)
{
  const char routineName[] = "scatterIntegral";
  register int  la, k, lap, kr, ip, mu, mup;

  char    filename[MAX_LINE_SIZE];
  bool_t  hunt, initialize, to_obs, to_obs_p;
  int     Np, Nread, Nwrite, ij, lamu;
  double *v_emit, v0, vN, *v_abs = NULL, *vp = NULL, *wv = NULL,
         *rii = NULL, *adamp, *Jbar, cDop, *RIInorm, *I = NULL,
        **Imup, *Ik, *Pj, *gamma, **v_los, *phi_emit, wmup, *sv;
  Atom *atom;
  AtomicLine *line;
  AtomicContinuum *continuum;

  /* --- Calculate the angle-dependent scattering integral when
         angle-dependent scattering is requested. --   -------------- */

  atom = PRDline->atom;

  if (!PRDline->PRD) {
    sprintf(messageStr, "Line %d -> %d of %2s is not a PRD line",
	    PRDline->j, PRDline->i, atom->ID);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }
  getCPU(3, TIME_START, NULL);

  cDop = (NM_TO_M * PRDline->lambda0) / (4.0 * PI);

  /* --- When called for the first time open temporary file for
         storage of redistribution weights --            ------------ */

  initialize = FALSE;
  if (PRDline->fp_GII == NULL) {
    sprintf(filename,
	    (atom->ID[1] == ' ') ? "PRD_%.1s_%d-%d.dat" : "PRD_%s_%d-%d.dat",
	    atom->ID, PRDline->j, PRDline->i);

    /* --- First try if file exists and can be opened for reading - - */

    if ((PRDline->fp_GII = fopen(filename, "r")) != NULL) {
      sprintf(messageStr,
	      "Using file %s with existing redistribution weights", filename);
      Error(WARNING, routineName, messageStr);
    } else {
      initialize = TRUE;
      if ((PRDline->fp_GII = fopen(filename, "w+")) == NULL) {
	sprintf(messageStr, "Unable to open temporary file %s", filename);
	Error(ERROR_LEVEL_2, routineName, messageStr);
      }
    }
  }
  if (!initialize) rewind(PRDline->fp_GII);

  /* --- Temporary storage space --                    -------------- */

  Imup    = matrix_double(PRDline->Nlambda, atmos.Nspace);
  Ik      = (double *) malloc(PRDline->Nlambda * sizeof(double));
  v_emit  = (double *) malloc(atmos.Nspace * sizeof(double));
  v_abs   = (double *) malloc(PRDline->Nlambda * sizeof(double));
  adamp   = (double *) malloc(atmos.Nspace * sizeof(double));
  v_los   = matrix_double(atmos.Nrays, atmos.Nspace);
  gamma   = (double *) malloc(atmos.Nspace * sizeof(double));
  Pj      = (double *) malloc(atmos.Nspace * sizeof(double));
  Jbar    = (double *) malloc(atmos.Nspace * sizeof(double));
  RIInorm = (double *) malloc(atmos.Nspace * sizeof(double));
  sv      = (double *) malloc(atmos.Nspace * sizeof(double));

  /* --- Evaluate first the total rate Pj out of the line's upper level
         and then the coherency fraction gamma --      -------------- */

  for (k = 0;  k < atmos.Nspace;  k++)
    Pj[k] =  PRDline->Qelast[k];
  for (ip = 0;  ip < atom->Nlevel;  ip++) {
    ij = ip * atom->Nlevel + PRDline->j;
    for (k = 0;  k < atmos.Nspace;  k++)
      Pj[k] += atom->C[ij][k];
  }
  for (kr = 0;  kr < atom->Nline;  kr++) {
    line = &atom->line[kr];
    if (line->j == PRDline->j)
      for (k = 0;  k < atmos.Nspace;  k++) Pj[k] += line->Rji[k];
    if (line->i == PRDline->j)
      for (k = 0;  k < atmos.Nspace;  k++) Pj[k] += line->Rij[k];
  }
  for (kr = 0;  kr < atom->Ncont;  kr++) {
    continuum = &atom->continuum[kr];
    if (continuum->j == PRDline->j)
      for (k = 0;  k < atmos.Nspace;  k++) Pj[k] += continuum->Rji[k];
    if (continuum->i == PRDline->j)
      for (k = 0;  k < atmos.Nspace;  k++) Pj[k] += continuum->Rij[k];
  }
  for (k = 0;  k < atmos.Nspace;  k++) {
    gamma[k] = atom->n[PRDline->i][k] / atom->n[PRDline->j][k] *
      PRDline->Bij / Pj[k];
  }
  /* --- Store line-of-sight velocity in Doppler units to avoid
         having to recompute it for every wavelength -- ------------- */

  for (mu = 0;  mu < atmos.Nrays;  mu++) {
    for (k = 0;  k < atmos.Nspace;  k++)
      v_los[mu][k] = vproject(k, mu) / atom->vbroad[k];
  }
  /* --- Depth-dependent damping parameter --          -------------- */

  for (k = 0;  k < atmos.Nspace;  k++) {
    Jbar[k]  = PRDline->Rij[k] / PRDline->Bij;
    sv[k]    = 1.0 / (SQRTPI * atom->vbroad[k]);
    adamp[k] =
      (PRDline->Grad + PRDline->Qelast[k]) * cDop / atom->vbroad[k];
  }
  /* --- Outer loop over emission wavelength --        -------------- */

  for (la = 0;  la < PRDline->Nlambda;  la++) {

    /* --- Loop over emission angle (down and up) --   -------------- */

    for (mu = 0;  mu < atmos.Nrays;  mu++) {
      for (to_obs = 0;  to_obs <= 1;  to_obs++) {
	lamu = 2*(atmos.Nrays*la + mu) + to_obs;

	if (atmos.moving ||
	    (PRDline->polarizable && input.StokesMode > FIELD_FREE))
	  phi_emit = PRDline->phi[lamu];
	else
	  phi_emit = PRDline->phi[la];

	for (k = 0;  k < atmos.Nspace;  k++) {
	  PRDline->rho_prd[lamu][k] = 0.0;
	  RIInorm[k] = 0.0;
	}
        /* --- Depth-dependent emission wavelengths for this
               direction (in Doppler units). --        -------------- */

	if (to_obs) {
	  for (k = 0;  k < atmos.Nspace;  k++)
	    v_emit[k] = (PRDline->lambda[la] - PRDline->lambda0) * CLIGHT /
	      (atom->vbroad[k] * PRDline->lambda0) + v_los[mu][k];
	} else {
	  for (k = 0;  k < atmos.Nspace;  k++)
	    v_emit[k] = (PRDline->lambda[la] - PRDline->lambda0) * CLIGHT /
	      (atom->vbroad[k] * PRDline->lambda0) - v_los[mu][k];
	}

        /* --- Loop over absorption directions --      -------------- */

	for (mup = 0;  mup < atmos.Nrays;  mup++) {
	  wmup = 0.5 * atmos.wmu[mup];
	  for (to_obs_p = 0;  to_obs_p <= 1;  to_obs_p++) {

	    /* --- Read specific intensity in this direction for all
                   wavelengths --                      -------------- */

            for (lap = 0;  lap < PRDline->Nlambda;  lap++)
	      readImu(PRDline->Nblue + lap, mup, to_obs_p, Imup[lap]);

	    /* --- Loop over space --                  -------------- */

	    for (k = 0;  k < atmos.Nspace;  k++) {
	      for (lap = 0;  lap < PRDline->Nlambda;  lap++) {
		Ik[lap] = Imup[lap][k];

		/* --- Array of absorption wavelengths in this
		       direction --                    -------------- */

		if (to_obs_p)
		  v_abs[lap] =
		    (PRDline->lambda[lap] - PRDline->lambda0) * CLIGHT /
		    (PRDline->lambda0 * atom->vbroad[k]) + v_los[mup][k];
		else
		  v_abs[lap] =
		    (PRDline->lambda[lap] - PRDline->lambda0) * CLIGHT /
		    (PRDline->lambda0 * atom->vbroad[k]) - v_los[mup][k];
	      }
              /* --- Setup spline interpolation coefficients for the
		     integration over absorption intensity -- ------- */

	      switch (representation) {
	      case SPLINE:
		splineCoef(PRDline->Nlambda, v_abs, Ik);
		break;
	      case EXP_SPLINE:
		exp_splineCoef(PRDline->Nlambda, v_abs, Ik, TENSION);
		break;
	      case LINEAR: break;
	      }
	      /* --- Establish integration limits over absorption
                     wavelength, using only regions where the
                     redistribution function is non-zero. -- -------- */
 
	      if (fabs(v_emit[k]) < PRD_QCORE) {
		v0 = -PRD_QWING;
		vN =  PRD_QWING;
	      } else {
		if (fabs(v_emit[k]) < PRD_QWING) {
		  if (v_emit[k] > 0.0) {
		    v0 = -PRD_QWING;
		    vN = v_emit[k] + PRD_QSPREAD;
		  } else {
		    v0 = v_emit[k] - PRD_QSPREAD;
		    vN = PRD_QWING;
		  }
		} else {
		  v0 = v_emit[k] - PRD_QSPREAD;
		  vN = v_emit[k] + PRD_QSPREAD;
		}
	      }
	      Np = (int) ((vN - v0) / PRD_DQ) + 1;
	      vp = (double *) realloc(vp,  Np * sizeof(double));
	      for (lap = 1, vp[0] = v0;  lap < Np;  lap++)
		vp[lap] = vp[lap - 1] + PRD_DQ;

	      /* --- Interpolate specific intensity onto fine grid.
                     Choose linear, spline or exponential spline
                     interpolation --                  -------------- */

	      I = (double *) realloc(I,  Np * sizeof(double));
	      switch (representation) {
	      case LINEAR:
		Linear(PRDline->Nlambda, v_abs, Ik, Np, vp, I, hunt=TRUE);
		break;
	      case SPLINE:
		splineEval(Np, vp, I, hunt=TRUE);
		break;
	      case EXP_SPLINE:
		exp_splineEval(Np, vp, I, hunt=TRUE);
		break;
	      }
	      /* --- Compute the redistribution weights -- ---------- */
	      
	      rii = (double *) realloc(rii, Np * sizeof(double));
	      if (initialize) {

		/* --- Integration weights (See: Press et al.
		       Numerical Recipes, p. 107, eq. 4.1.12) -- ---- */

		wv = (double *) realloc(wv,  Np * sizeof(double));
		wv[0] = 5.0/12.0  * PRD_DQ;
		wv[1] = 13.0/12.0 * PRD_DQ;
		for (lap = 2;  lap < Np-2;  lap++) wv[lap] = PRD_DQ;
		wv[Np-1] = 5.0/12.0  * PRD_DQ;
		wv[Np-2] = 13.0/12.0 * PRD_DQ;

		/* --- RII(x, mu, x', mu')/phi(x, mu) -- ------------ */

		for (lap = 0;  lap < Np;  lap++) {
		  rii[lap] = RII(v_emit[k], vp[lap], adamp[k], mu, mup) *
		    (sv[k] / phi_emit[k]) * wv[lap] * wmup;
		}
		if ((Nwrite = fwrite(rii, sizeof(double), Np,
				     PRDline->fp_GII)) != Np) {
		  sprintf(messageStr,
		"Unable to write proper number of redistribution weights\n"
		" Wrote %d instead of %d.\n Line %d -> %d, la = %d, k = %d",
			  Nwrite, Np, PRDline->j, PRDline->i, la, k);
		  Error(ERROR_LEVEL_2, routineName, messageStr);
		}
	      } else {
		if ((Nread = fread(rii, sizeof(double), Np,
				   PRDline->fp_GII)) != Np) {
		  sprintf(messageStr,
		"Unable to read proper number of redistribution weights\n"
	        " Read %d instead of %d.\n Line %d -> %d, la = %d, k = %d",
			  Nread, Np, PRDline->j, PRDline->i, la, k);
		  Error(ERROR_LEVEL_2, routineName, messageStr);
		}
	      }
	      /* --- Inner wavelength loop doing actual wavelength
                     integration over absorption wavelengths -- ----- */

	      for (lap = 0;  lap < Np;  lap++) {
		RIInorm[k]                += rii[lap];
		PRDline->rho_prd[lamu][k] += I[lap] * rii[lap];
	      }
	    }
	  }
	}

	for (k = 0;  k < atmos.Nspace;  k++) {
	  PRDline->rho_prd[lamu][k] = 1.0 +
	    gamma[k] * (PRDline->rho_prd[lamu][k]/RIInorm[k] - Jbar[k]);
	}
      }
    }
  }
  /* --- Clean temporary variable space --             -------------- */

  freeMatrix((void **) Imup);
  freeMatrix((void **) v_los);

  free(Ik);     free(v_abs);  free(v_emit);  free(adamp);
  free(gamma);  free(Pj);
  free(wv);     free(I);       free(rii);
  free(vp);     free(Jbar);   free(RIInorm); free(sv);

  sprintf(messageStr, "Scatter Int %5.1f", PRDline->lambda0);
  getCPU(3, TIME_POLL, messageStr);
}
/* ------- end ---------------------------- PRDAngleScatter.c ------- */
