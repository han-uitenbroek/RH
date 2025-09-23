/* ------- file: -------------------------- formal.c ----------------

       Version:       rh2.0, 2-D Cartesian
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Tue Sep 23 10:50:09 2025 --

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
#include "error.h"
#include "background.h"
#include "inputs.h"
#include "constant.h"
#include "statistics.h"
#include "xdr.h"

/* --- Function prototypes --                          -------------- */

void Feautrier_2D(int nspect, double *chi, double *S,
		  enum FeautrierOrder F_order, double *I, double *Psi);


/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern Geometry geometry;
extern Spectrum spectrum;
extern InputData input;
extern char   messageStr[];


/* ------- begin -------------------------- Formal.c ---------------- */

double Formal(int nspect, bool_t eval_operator, bool_t redistribute)
{
  const char routineName[] = "Formal";
  register int k, l, mu, n;

  bool_t   initialize, boundbound, polarized_as, polarized_c,
           PRD_angle_dep, to_obs, solveStokes, angle_dep;
  int      Nspace = atmos.Nspace, Nrays = atmos.Nrays, lamuk;
  long int idx, idx0;
  double  *I, *chi, *S, **Ipol, **Spol, *Psi, *Jdag, wmu, dJmax, dJ,
          *eta_Q, *eta_U, *eta_V, *eta_c_Q, *eta_c_U, *eta_c_V,
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
    Spol = matrix_double(4, Nspace);
    S    = Spol[0];
    Ipol = matrix_double(4, Nspace);
    I    = Ipol[0];
  } else {
    S = (double *) malloc(Nspace * sizeof(double));
    I = (double *) malloc(Nspace * sizeof(double));
  }
  chi  = (double *) malloc(Nspace * sizeof(double));

  /* --- Store current mean intensity, initialize new one to zero - - */

  Jdag = (double *) malloc(Nspace * sizeof(double));
  if (input.limit_memory) {
    J = (double *) malloc(atmos.Nspace * sizeof(double));
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
      wmu = 0.5 * geometry.wmu[mu];
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
	  
	  if (input.S_interpolation_stokes == DELO_PARABOLIC) {
	    Piece_Stokes_2D(&geometry, nspect, mu, to_obs,
			    chi, Spol, Ipol, Psi);
	  } else if (input.S_interpolation_stokes == DELO_BEZIER3) {
	    Piece_Stokes_Bezier3_2D(&geometry, nspect, mu, to_obs,
			    chi, Spol, Ipol, Psi);
	  } else {sprintf(messageStr,
			  "Unknown polarization solver: %d",
			  input.S_interpolation_stokes);
	    Error(ERROR_LEVEL_1, routineName, messageStr);
	  }
	} else {
	  for (k = 0;  k < Nspace;  k++)  S[k] /= chi[k];

	  /* --- Intensity only --                     -------------- */

	  if (input.S_interpolation == S_LINEAR) {
	    Piecewise_Linear_2D(&geometry, nspect, mu, to_obs,
				chi, S, I, Psi);
	  } else if (input.S_interpolation == S_PARABOLIC) {
	    Piecewise_2D(&geometry, nspect, mu, to_obs,
			 chi, S, I, Psi);
	  } else if (input.S_interpolation == S_BEZIER3) {
	    Piecewise_Bezier3_2D(&geometry, nspect, mu, to_obs,
				chi, S, I, Psi);
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
	  
	  /* --- Accumulate mean intensity --        ---------------- */

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

      for (l = 0;  l < geometry.Nx;  l++)
	spectrum.I[nspect*Nrays + mu][l] = I[l];

      if (solveStokes) {
	for (l = 0;  l < geometry.Nx;  l++) {
	  spectrum.Stokes_Q[nspect*Nrays + mu][l] = Ipol[1][l];
	  spectrum.Stokes_U[nspect*Nrays + mu][l] = Ipol[2][l];
	  spectrum.Stokes_V[nspect*Nrays + mu][l] = Ipol[3][l];
	}
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
    if (input.Eddington) {
      Feautrier_2D(nspect, chi, S, STANDARD, I, Psi);
      if (eval_operator) {
	for (k = 0;  k < Nspace;  k++) Psi[k] /= chi[k];
	addtoGamma(nspect, geometry.wmu[0], I, Psi);
      }
      if (spectrum.updateJ) {
	for (k = 0;  k < Nspace;  k++)
	  spectrum.J[nspect][k] += I[k] * geometry.wmu[0];
	addtoRates(nspect, 0, 0, geometry.wmu[0], I, redistribute);
      }
    } else {
      for (mu = 0;  mu < Nrays;  mu++) {
	wmu = 0.5 * geometry.wmu[mu];
	for (to_obs = 0;  to_obs <= 1;  to_obs++) {
	  if (input.S_interpolation == S_LINEAR) {
	    Piecewise_Linear_2D(&geometry, nspect, mu, to_obs,
				chi, S, I, Psi);
	  } else if (input.S_interpolation == S_PARABOLIC) {
	    Piecewise_2D(&geometry, nspect, mu, to_obs,
			 chi, S, I, Psi);
	  } else if (input.S_interpolation == S_BEZIER3) {
	    Piecewise_Bezier3_2D(&geometry, nspect, mu, to_obs,
				chi, S, I, Psi);
	  } else {
	    sprintf(messageStr,
		    "Unknown radiation solver: %d",
		    input.S_interpolation);
	    Error(ERROR_LEVEL_1, routineName, messageStr);
	  }
	  if (eval_operator) {
	    for (k = 0;  k < Nspace;  k++) Psi[k] /= chi[k];
	    addtoGamma(nspect, wmu, I, Psi);
	  }
          if (spectrum.updateJ) {
	    for (k = 0;  k < Nspace;  k++)
	      J[k] += I[k] * wmu;
	    addtoRates(nspect, mu, to_obs, wmu, I, redistribute);

	    /* --- Accumulate gas-frame mean intensity, which is the same
	           as J in the angle-independent case --   ------------ */
	
	    if (atmos.NPRDactive > 0  &&
		input.PRD_angle_dep == PRD_ANGLE_APPROX) {
	      for (k = 0;  k < Nspace;  k++)
		spectrum.Jgas[nspect][k] += I[k] * geometry.wmu[mu];
	    }
	  }
	}
	for (l = 0;  l < geometry.Nx;  l++)
	  spectrum.I[nspect*Nrays + mu][l] = I[l];
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

/* ------- begin -------------------------- Feautrier_2D.c ---------- */

#define  A_SIXTH  0.166666666667


void Feautrier_2D(int nspect, double *chi, double *S,
		  enum FeautrierOrder F_order, double *I, double *Psi)
{
  const char routineName[] = "Feautrier_2D";
  register int k, l, kl;

  int     Ndep = geometry.Nz;
  double  r0, h0, rN, hN, f0, fN, Ak, Ck, tau0 = 0.0, Bnu[2], zmu,
          dtau_mid, Iplus, *dtau, *abc, *A1, *C1, *F, *G, *Q,
         *Stmp, *ztmp, T[2], *chi_l, *S_l, *I_l, *Psi_l;

  dtau = (double *) malloc(geometry.Nz * sizeof(double));
  abc  = (double *) malloc(geometry.Nz * sizeof(double));
  Q    = (double *) malloc(geometry.Nz * sizeof(double));
  A1   = (double *) malloc(geometry.Nz * sizeof(double));
  C1   = (double *) malloc(geometry.Nz * sizeof(double));
  F    = (double *) malloc(geometry.Nz * sizeof(double));
  G    = (double *) malloc(geometry.Nz * sizeof(double));
  ztmp = (double *) malloc(geometry.Nz * sizeof(double));
  Stmp = (double *) malloc(geometry.Nz * sizeof(double));

  chi_l = (double *) malloc(geometry.Nz * sizeof(double));
  S_l   = (double *) malloc(geometry.Nz * sizeof(double));
  I_l   = (double *) malloc(geometry.Nz * sizeof(double));
  if (Psi) Psi_l = (double *) malloc(geometry.Nz * sizeof(double));

  zmu = 0.5 / geometry.muz[0];

  for (l = 0;  l < geometry.Nx;  l++) {
    for (k = 0, kl = l;  k < geometry.Nz;  k++, kl += geometry.Nx) {
      chi_l[k] = chi[kl];
      S_l[k] = S[kl];
    }
    for (k = 0;  k < geometry.Nz-1;  k++)
      dtau[k] = zmu * (chi_l[k] + chi_l[k+1]) * geometry.dz[k];


    /* --- Upper boundary condition:  I^- = r0*I^+ + h0 -- ---------- */

    switch (geometry.bvalue[TOP]) {
    case IRRADIATED:
      r0 = 0.0;
      h0 = geometry.Itop[nspect][l];
      break;
    case ZERO:
      r0 = h0 = 0.0;
      break;
    case THERMALIZED:
      r0 = 0.0;
      T[0] = atmos.T[l];
      T[1] = atmos.T[geometry.Nx + l];

      Planck(2, T, spectrum.lambda[nspect], Bnu);
      h0 = Bnu[0] - (Bnu[1] - Bnu[0]) / dtau[0];
      break;
    case REFLECTIVE:
      Error(ERROR_LEVEL_2, routineName, 
	    "Boundary condition REFLECTIVE not implemented");
      break;
    }

    f0      = (1.0 - r0) / (1.0 + r0);
    abc[0]  = 1.0 + 2.0*f0 / dtau[0];
    C1[0]   = 2.0 / SQ(dtau[0]);
    Stmp[0] = S_l[0] + 2.0*h0 / ((1.0 + r0)*dtau[0]);
    if (F_order == FEAUTRIER_HERMITE) {
      C1[0]   -= 2.0*A_SIXTH;
      Stmp[0] += 2.0*A_SIXTH * (S_l[1] - S_l[0]);
    }
    /* --- Lower boundary condition:  I^+ = rN*I^- + hN -- ---------- */

    switch (geometry.bvalue[BOTTOM]) {
    case IRRADIATED:
      rN = 0.0;
      hN = geometry.Ibottom[nspect][l];
      break;
    case ZERO:
      rN = hN = 0.0;
      break;
    case THERMALIZED:
      rN = 0.0;
      T[0] = atmos.T[(geometry.Nz-2)*geometry.Nx + l];
      T[1] = atmos.T[(geometry.Nz-1)*geometry.Nx + l];

      Planck(2, T, spectrum.lambda[nspect], Bnu);
      hN = Bnu[1] - (Bnu[0] - Bnu[1]) / dtau[Ndep-2];
      break;
    case REFLECTIVE:
      Error(ERROR_LEVEL_2, routineName, 
	    "Boundary condition REFLECTIVE not implemented");
      break;
    }

    fN           = (1.0 - rN) / (1.0 + rN);
    abc[Ndep-1]  = 1.0 + 2.0*fN / dtau[Ndep-2];
    A1[Ndep-1]   = 2.0 / SQ(dtau[Ndep-2]);
    Stmp[Ndep-1] = S_l[Ndep-1] + 2.0*hN / ((1.0 + rN)*dtau[Ndep-2]);
    if (F_order == FEAUTRIER_HERMITE) {
      A1[Ndep-1]   -= 2.0*A_SIXTH;
      Stmp[Ndep-1] += 2.0*A_SIXTH * (S_l[Ndep-2] - S_l[Ndep-1]);
    }

    for (k = 1;  k < Ndep-1;  k++) {
      dtau_mid = 0.5*(dtau[k] + dtau[k-1]);
      A1[k]   = 1.0 / (dtau_mid * dtau[k-1]);
      C1[k]   = 1.0 / (dtau_mid * dtau[k]);
      abc[k]  = 1.0;
      Stmp[k] = S_l[k];
    }
    if (F_order == FEAUTRIER_HERMITE) {
      for (k = 1;  k < Ndep-1;  k++) {
	Ak     = A_SIXTH * (1.0 - 0.5*SQ(dtau[k])*A1[k]);
	Ck     = A_SIXTH * (1.0 - 0.5*SQ(dtau[k-1])*C1[k]);
	A1[k] -= Ak;
	C1[k] -= Ck;
	Stmp[k] += Ak*(S_l[k-1] - S_l[k]) + Ck*(S_l[k+1] - S_l[k]);
      }
    }
    /* --- Start the elimination --                    -------------- */

    F[0]    = abc[0] / C1[0];
    ztmp[0] = Stmp[0] / (abc[0] + C1[0]);
    for (k = 1;  k < Ndep-1;  k++) {
      F[k]    = (abc[k] + A1[k]*F[k-1]/(1.0 + F[k-1])) / C1[k];
      ztmp[k] = (Stmp[k] + A1[k]*ztmp[k-1]) / (C1[k] * (1.0 + F[k]));
    }
    /* --- Now backsubstitution --                     --- ---------- */

    I_l[Ndep-1] = (Stmp[Ndep-1]+ A1[Ndep-1]*ztmp[Ndep-2]) /
      (abc[Ndep-1] + A1[Ndep-1]*(F[Ndep-2] / (1.0 + F[Ndep-2])));
    for (k = Ndep-2;  k >= 0;  k--)
      I_l[k] = I_l[k+1] / (1.0 + F[k]) + ztmp[k];

    /* --- If necessary evaluate the diagonal operator -- ----------- */

    if (Psi) {
      if (F_order == FEAUTRIER_HERMITE) {
	sprintf(messageStr,
	"Higher order for diagonal operator calculation not yet implemented");
	Error(ERROR_LEVEL_1, routineName, messageStr);
      }

      G[Ndep-1] = abc[Ndep-1] / A1[Ndep-1];
      for (k = Ndep-2;  k >= 1;  k--)
	G[k] = (abc[k] + C1[k]*G[k+1]/(1.0 + G[k+1])) / A1[k];

      Psi_l[0] = 1.0 / (abc[0] + C1[0]*G[1]/(1.0 + G[1]));
      for (k = 1;  k < Ndep-1;  k++)
	Psi_l[k] = 1.0 / (abc[k] + A1[k]*F[k-1]/(1.0 + F[k-1]) +
			  C1[k]*G[k+1]/(1.0 + G[k+1]));
      Psi_l[Ndep-1] =
	1.0 / (abc[Ndep-1] + A1[Ndep-1]*F[Ndep-2]/(1.0 + F[Ndep-2]));
    }
    /* --- Emergent intensity --                       -------------- */

    Iplus = (1.0 + f0)*I_l[0] - h0/(1.0 + r0);
    if (tau0) 
      spectrum.I[nspect][l] = (Iplus - S_l[0])*exp(-tau0) + S_l[0];
    else
      spectrum.I[nspect][l] = Iplus;

    for (k = 0, kl = l;  k < geometry.Nz;  k++, kl += geometry.Nx)
      I[kl] = I_l[k];
    if (Psi) {
      for (k = 0, kl = l;  k < geometry.Nz;  k++, kl += geometry.Nx)
	Psi[kl] = Psi_l[k];
    }
  }

  free(dtau);  free(abc);  free(Q);   free(A1);
  free(C1);    free(F);    free(G);   free(ztmp);
  free(Stmp);

  free(chi_l);  free(S_l);  free(I_l);
  if (Psi) free(Psi_l);
}
/* ------- end ---------------------------- Feautrier_2D.c ---------- */

/* ------- begin -------------------------- Hydrostatic.c ----------- */

void Hydrostatic(int NmaxIter, double iterLimit)
{
  const char routineName[] = "Hydrostatic";

  if (atmos.hydrostatic) {
    sprintf(messageStr,
	    "Can only establish hydrostatic equilibrium in 1-D geometry");
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }
}
/* ------- end ---------------------------- Hydrostatic.c ----------- */
