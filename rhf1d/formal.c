/* ------- file: -------------------------- formal_1D.c --------------

       Version:       rh1.0, 1-D plane-parallel
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Tue Feb 20 15:18:13 2024 --

       --------------------------                      ----------RH-- */

/* --- Formal solution with given source function, and allowing for
       a PRD emission profile, polarized line radiation, polarized
       background lines, and scattering background polarization -- -- */

 
#include <stdlib.h>
#include <math.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "geometry.h"
#include "spectrum.h"
#include "constant.h"
#include "background.h"
#include "inputs.h"
#include "error.h"
#include "xdr.h"
#include "bezier.h"


/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern Geometry geometry;
extern Spectrum spectrum;
extern InputData input;
extern char messageStr[];


/* ------- begin -------------------------- Formal.c ---------------- */

double Formal(int nspect, bool_t eval_operator, bool_t redistribute)
{
  const char routineName[] = "Formal";
  register int k, mu, n;

  bool_t   initialize, boundbound, polarized_as, polarized_c,
           PRD_angle_dep, to_obs, solveStokes, angle_dep;
  enum     FeautrierOrder F_order;     
  int      Nspace = atmos.Nspace, Nrays = atmos.Nrays, lamuk;
  long int idx,idx0;
  double  *I, *chi, *S, **Ipol, **Spol, *Psi, *Jdag, wmu, dJmax, dJ,
          *J20dag, musq, threemu1, threemu2, *J, *J20, *lambda, sign,
           lambda_gas, lambda_prv, lambda_nxt, fac, dl, frac;
  ActiveSet *as;

  /* --- Retrieve active set as of transitions at wavelength nspect - */

  as = &spectrum.as[nspect];
  alloc_as(nspect, eval_operator);
  
  /* --- Check whether current active set includes a bound-bound
         and/or polarized transition and/or angle-dependent PRD
         transition, and/or polarization through background scattering.
         Otherwise, only angle-independent opacity and source functions
         are needed --                                 -------------- */ 

  /* --- Check for bound-bound transition in active set -- ---------- */

  boundbound    = containsBoundBound(as);

  /* --- Check for line with angle-dependent PRD in set -- ---------- */

  PRD_angle_dep = (containsPRDline(as) &&
		   (input.PRD_angle_dep == PRD_ANGLE_DEP  ||
		    input.PRD_angle_dep == PRD_ANGLE_APPROX));

  /* --- Check for polarized bound-bound transition in active set - - */

  polarized_as  = containsPolarized(as);

  /* --- Check for polarized bound-bound transition in background - - */

  polarized_c   = atmos.backgrflags[nspect].ispolarized;

  /* --- Determine if we solve for I, or for I, Q, U, V -- ---------- */

  solveStokes   = (input.StokesMode == FULL_STOKES &&
		   (polarized_as || polarized_c || input.backgr_pol));

  /* --- Determine if we have to do angle-dependent opacity and
         emissivity --                                 -------------- */

  angle_dep     = (polarized_as || polarized_c
		   || PRD_angle_dep || input.backgr_pol ||
		   (atmos.moving &&
		    (boundbound || atmos.backgrflags[nspect].hasline)));

  /* --- Allocate temporary space --                   -------------- */

  if (eval_operator)
    Psi = (double *) malloc(Nspace * sizeof(double));
  else
    Psi = NULL;

  if (solveStokes) {
    Ipol = matrix_double(4, Nspace);
    I    = Ipol[0];
    Spol = matrix_double(4, Nspace);
    S    = Spol[0];
  } else {
    I = (double *) malloc(Nspace * sizeof(double));
    S = (double *) malloc(Nspace * sizeof(double));
  }
  chi = (double *) malloc(Nspace * sizeof(double));

  /* --- Store current mean intensity, initialize new one to zero - - */

  Jdag = (double *) malloc(Nspace * sizeof(double));
  if (input.limit_memory) {
    J = (double *) malloc(atmos.Nspace *sizeof(double));
    readJlambda(nspect, Jdag);
  } else {
    J = spectrum.J[nspect];
    for (k = 0;  k < Nspace;  k++) Jdag[k] = J[k];
  }
  if (spectrum.updateJ)
    for (k = 0;  k < Nspace;  k++) J[k] = 0.0;

  /* --- Store current anisotropy, initialize new one to zero ---- -- */

  if (input.backgr_pol) {
    J20dag = (double *) malloc(Nspace * sizeof(double));
    if (input.limit_memory) {
      J20 = (double *) malloc(Nspace * sizeof(double));
      readJ20lambda(nspect, J20dag);
    } else {
      J20 = spectrum.J20[nspect];
      for (k = 0;  k < Nspace;  k++)
	J20dag[k] = J20[k];
    }
    if (spectrum.updateJ)
      for (k = 0;  k < Nspace;  k++) J20[k] = 0.0;
  }
  /* --- Case of angle-dependent opacity and source function -- ----- */

  if (angle_dep) {
    for (mu = 0;  mu < Nrays;  mu++) {
      wmu  = 0.5 * geometry.wmu[mu];
      if (input.backgr_pol) {
	musq = SQ(geometry.muz[mu]);
	threemu1 = TWOSQRTTWO * (3.0*musq - 1.0);
	threemu2 = (3.0 * TWOSQRTTWO) * (musq - 1.0);
      }
      for (to_obs = 0;  to_obs <= 1;  to_obs++) {
	initialize = (mu == 0 && to_obs == 0);

	if (initialize || atmos.backgrflags[nspect].hasline)
	  readBackground(nspect, mu, to_obs);

	if (initialize || boundbound)
	  Opacity(nspect, mu, to_obs, initialize);

	if (eval_operator) addtoCoupling(nspect);
	for (k = 0;  k < Nspace;  k++) {
	  chi[k] = as->chi[k] + as->chi_c[k];
	  S[k]   = as->eta[k] + as->eta_c[k] + as->sca_c[k]*Jdag[k];
	}
	if (input.backgr_pol) {
	  for (k = 0;  k < Nspace;  k++) {
	    S[k] += threemu1 * as->sca_c[k]*J20dag[k];
	  }
	} 
	if (solveStokes) {
	  for (k = Nspace;  k < 4*Nspace;  k++) Spol[0][k] = 0.0;

          /* --- Add emissivity due to active set for Q, U, V -- ---- */

          if (polarized_as) {
            for (k = Nspace;  k < 4*Nspace;  k++)
	      Spol[0][k] += as->eta[k];
	  }
          /* --- Add emissivity due to background lines -- ---------- */

          if (polarized_c) {
            for (k = Nspace;  k < 4*Nspace;  k++)
	      Spol[0][k] += as->eta_c[k];
	  }
	  /* --- Add emissivity due to background scattering -- ----- */

          if (input.backgr_pol) {
            for (k = 0;  k < Nspace;  k++) {
              Spol[1][k] += threemu2 * as->sca_c[k]*J20dag[k];
	    }
	  }
	  for (n = 0;  n < 4;  n++) {
	    for (k = 0;  k < Nspace;  k++)
	      Spol[n][k] /= chi[k];
	  }
	  /* --- Polarized transfer --                 -------------- */
	  
	  if (input.S_interpolation_stokes == DELO_BEZIER3) {
	    Piece_Stokes_Bezier3_1D(nspect, mu, to_obs,
				    chi, Spol, Ipol, Psi);
	  } else if (input.S_interpolation_stokes == DELO_PARABOLIC) {
	    Piece_Stokes_1D(nspect, mu, to_obs, chi, Spol, Ipol, Psi);
	  } else { sprintf(messageStr,
		    "Unknown polarization solver: %d",
		    input.S_interpolation_stokes);
	    Error(ERROR_LEVEL_1, routineName, messageStr);	    
	  }
	} else {
	  
	  for (k = 0;  k < Nspace;  k++)
	    S[k] /= chi[k];

	  /* --- Intensity only --                     -------------- */
	  
	  if (input.S_interpolation == S_LINEAR) {
	    Piecewise_Linear_1D(nspect, mu, to_obs, chi, S, I, Psi);
	  } else if (input.S_interpolation == S_PARABOLIC) {
	    Piecewise_1D(nspect, mu, to_obs, chi, S, I, Psi);
	  } else if (input.S_interpolation == S_BEZIER3) {
	    Piecewise_Bezier3_1D(nspect, mu, to_obs, chi, S, I, Psi);
	  } else {
	    sprintf(messageStr,
		    "Unknown radiation solver: %d",
		    input.S_interpolation);
	    Error(ERROR_LEVEL_1, routineName, messageStr);
	  }
	}
	
	if (eval_operator) {
          for (k = 0;  k < Nspace;  k++) Psi[k] /= chi[k];
	  addtoGamma(nspect, wmu, I, Psi);
	}

        if (spectrum.updateJ) {

	  /* --- Accumulate mean intensity and rates -- ----------- */

	  for (k = 0;  k < Nspace;  k++)
	    J[k] += wmu * I[k];
	  addtoRates(nspect, mu, to_obs, wmu, I, redistribute);

	  /* --- Accumulate anisotropy --            -------------- */

	  if (input.backgr_pol) {
            if (solveStokes) {
	      for (k = 0;  k < Nspace;  k++) {
		J20[k] += (threemu1 * Ipol[0][k] +
			   threemu2 * Ipol[1][k]) * wmu;
	      }
	    } else {
	      for (k = 0;  k < Nspace;  k++) {
		J20[k] += threemu1 * I[k] * wmu;
	      }
	    }
	  }

	  /* --- Accumulate gas-frame mean intensity --   ----------- */

	  if (atmos.NPRDactive > 0  &&
	      input.PRD_angle_dep == PRD_ANGLE_APPROX  &&
	      atmos.Nrays > 1) {
	    
	    if (input.prdh_limit_mem) {
	      sign = (to_obs) ? 1.0 : -1.0;

	      for (k = 0;  k < Nspace;  k++) {

		/* --- Observer's frame wavelenght grid -- ---------- */
		
		lambda = spectrum.lambda;

		/* -- Previous, current and next wavelength shifted to
		      gas rest frame --                -------------- */
		
		fac = (1.0 + spectrum.v_los[mu][k] * sign/CLIGHT);
		lambda_prv = lambda[ MAX(nspect-1,0)                 ] * fac;
		lambda_gas = lambda[ nspect                          ] * fac;
		lambda_nxt = lambda[ MIN(nspect+1,spectrum.Nspect-1) ] * fac;

		/* ---  Do lambda_prv and lambda_gas bracket
		        lambda points? --              -------------- */
		
		if (lambda_prv != lambda_gas) {

		  dl = lambda_gas - lambda_prv;
		  for (idx = 0;  idx < spectrum.Nspect;  idx++) {
		    if (lambda[idx] > lambda_prv  &&
			lambda[idx] <= lambda_gas) {
		      frac = (lambda[idx]-lambda_prv) / dl;
		      spectrum.Jgas[idx][k] += frac * wmu * I[k];
		    }
		  }
		} else {

		  /* --- Edge case, use constant extrapolation for
		         lambda[idx] < lambda gas --   -------------- */
		  
		  for (idx = 0;  idx < spectrum.Nspect;  idx++) {
		    if (lambda[idx] < lambda_gas)
		      spectrum.Jgas[idx][k] += wmu * I[k];
		  }
		}

		/* --- Do lambda_gas and lambda_nxt bracket
		       lambda points? --               -------------- */
		
		if (lambda_gas != lambda_nxt) {
		  dl = lambda_nxt - lambda_gas;
		  
		  for (idx = 0;  idx < spectrum.Nspect;  idx++) {
		    if (lambda[idx] > lambda_gas &&
			lambda[idx] < lambda_nxt) {
		      frac = (lambda[idx] - lambda_gas) / dl;
		      spectrum.Jgas[idx][k] += (1.0 - frac) * wmu * I[k];
		    }
		  }
		} else {
		  
		  /* --- Edge case, use constant extrapolation for
		         lambda[idx] > lambda gas --   -------------- */
		  
		  for (idx = 0;  idx < spectrum.Nspect;  idx++) {
		    if (lambda[idx] >  lambda_gas)
		      spectrum.Jgas[idx][k] += wmu * I[k];
		  }
		}
	      }
	    } else {

	      for (k = 0;  k < Nspace;  k++)  {
		lamuk = nspect * (atmos.Nrays * 2 * Nspace)
		  + mu * (2 * Nspace) + to_obs * (Nspace) + k;

		idx0 = (lamuk==0) ? 0 : spectrum.nc[lamuk-1];

		for ( idx = idx0;  idx <  spectrum.nc[lamuk];  idx++ )
		  spectrum.Jgas[spectrum.iprdh[idx]][k] +=
		    wmu *  spectrum.cprdh[idx] * I[k];
	      }
	    }
	  }
	  
	  if (containsPRDline(as) &&
	      input.PRD_angle_dep == PRD_ANGLE_DEP)
	    writeImu(nspect, mu, to_obs, I);
	}
      }
      /* --- Save emergent intensity --              -------------- */

      spectrum.I[nspect][mu] = I[0];
      if (solveStokes) {
	spectrum.Stokes_Q[nspect][mu] = Ipol[1][0];
	spectrum.Stokes_U[nspect][mu] = Ipol[2][0];
	spectrum.Stokes_V[nspect][mu] = Ipol[3][0];
      }
    }
  } else {

    /* --- The angle-independent case --               -------------- */

    readBackground(nspect, 0, 0);
    Opacity(nspect, 0, 0, initialize=TRUE);
    if (eval_operator) addtoCoupling(nspect);

    for (k = 0;  k < Nspace;  k++) {
      chi[k] = as->chi[k] + as->chi_c[k];
      S[k]   = (as->eta[k] +
		as->eta_c[k] + as->sca_c[k]*Jdag[k]) / chi[k];
    }
    for (mu = 0;  mu < Nrays;  mu++) {
      spectrum.I[nspect][mu] = Feautrier(nspect, mu, chi, S,
					 F_order=STANDARD, I, Psi);
      if (eval_operator) {
	for (k = 0;  k < Nspace;  k++) Psi[k] /= chi[k];
	addtoGamma(nspect, geometry.wmu[mu], I, Psi);
      }

      if (spectrum.updateJ) {
	for (k = 0;  k < Nspace;  k++) J[k] += I[k] * geometry.wmu[mu];
	addtoRates(nspect, mu, 0, geometry.wmu[mu], I, redistribute);

	/* --- Accumulate gas-frame mean intensity, which is the same
	       as J in the angle-independent case --   -------------- */
	
	if (atmos.NPRDactive > 0  &&
	    input.PRD_angle_dep == PRD_ANGLE_APPROX) {
	  for (k = 0;  k < Nspace;  k++)
	    spectrum.Jgas[nspect][k] += I[k] * geometry.wmu[mu];
	}
      }
    }
  }
  /* --- Write new J for current position in the spectrum -- -------- */

  dJmax = 0.0;
  if (spectrum.updateJ) {
    for (k = 0;  k < Nspace;  k++) {
      dJ = fabs(1.0 - Jdag[k]/J[k]);
      dJmax = MAX(dJmax, dJ);
    }
    if (input.limit_memory) {
      writeJlambda(nspect, J);
      if (input.backgr_pol)
	writeJ20lambda(nspect, J20);
    }
  }
  /* --- Clean up --                                 ---------------- */

  free_as(nspect, eval_operator);
  if (eval_operator) free(Psi);

  free(chi); 
  if (solveStokes) {
    freeMatrix((void **) Ipol);
    freeMatrix((void **) Spol);
  } else {
    free(I);
    free(S);
  }

  free(Jdag);
  if (input.limit_memory) free(J);
  if (input.backgr_pol) {
    free(J20dag);
    if (input.limit_memory) free(J20);
  }

  return dJmax;
}   
/* ------- end ---------------------------- Formal.c ---------------- */
