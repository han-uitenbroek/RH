/* ------- file: -------------------------- initial_xdr.c -----------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Mon Sep 22 15:54:37 2025 --

       --------------------------                      ----------RH-- */

/* --- Reads and/or computes the initial solution (populations and/or
       mean intensity J).

       XDR (external data representation) version.

       Possible options:

         LTE_POPULATIONS    -- Assume LTE populations initially.
         ZERO_RADIATION     -- Solve statistical equilibrium with
                               zero radiation field
         OLD_POPULATIONS    -- Read old populations from file
         OLD_POPS_AND_J     -- Read both old populations and J from file
         ESCAPE_PROBABILITY -- Not yet implemented
         OLD_J              -- Use mean intensities from previous solution
                               (Only implemented for wavelength_table).

       --                                              -------------- */


#include <fcntl.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "spectrum.h"
#include "accelerate.h"
#include "constant.h"
#include "statistics.h"
#include "error.h"
#include "inputs.h"
#include "xdr.h"


/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern Spectrum spectrum;
extern InputData input;
extern CommandLine commandline;
extern char messageStr[];

extern enum Topology topology;


/* ------- begin -------------------------- initSolution.c ---------- */

void initSolution(void)
{
  const char routineName[] = "initSolution";
  register int k, i, ij, nspect, n, kr, nact, mu;

  char     permission[3];
  bool_t   result, openJfile;
  int      la, j, niter, Nsr, Nplane, index, status, oflag,
           to_obs, lamuk, sign, ncoef, ilow, Nlamu, lamu;
  double   gijk, wla, twohnu3_c2, hc_k, twoc, fourPI, *J, *J20,
           *lambda, fac, lambda_prv, lambda_gas, lambda_nxt, dl,
           frac, lag;
  long int idx, lc;
  Atom     *atom;
  Molecule *molecule;
  ActiveSet *as;
  AtomicLine *line;
  AtomicContinuum *continuum;
  XDR xdrs;
  FILE *fp;

  getCPU(2, TIME_START, NULL);

  /* --- Allocate space for angle-averaged mean intensity -- -------- */

  if (!input.limit_memory) {
    spectrum.J = matrix_double(spectrum.Nspect, atmos.Nspace);

    /* --- If we do background polarization we need space for the
           anisotropy --                               -------------- */

    if (input.backgr_pol)
      spectrum.J20 = matrix_double(spectrum.Nspect, atmos.Nspace);
  }

  /* --- For the PRD angle approximation we need to store J in
         the gas frame --                              -------------- */
  
  if (input.PRD_angle_dep == PRD_ANGLE_APPROX  &&
      atmos.NPRDactive > 0) {

    spectrum.Jgas  = matrix_double(spectrum.Nspect, atmos.Nspace);
    spectrum.v_los = matrix_double(atmos.Nrays, atmos.Nspace);

    /* --- Calculate line of sight velocity --         -------------- */
    
    for (mu = 0;  mu < atmos.Nrays;  mu++) {
      for (k = 0;  k < atmos.Nspace;  k++) {
	spectrum.v_los[mu][k] = vproject(k, mu); // / vbroad[k];
      }
    }

    /* --- Precompute prd_rho interpolation coefficients if requested */
    
    if (!input.prdh_limit_mem) {

      for (nact = 0;  nact < atmos.Nactiveatom;  nact++) {
	atom = atmos.activeatoms[nact];

	for (kr = 0;  kr < atom->Nline;  kr++) {
	  line = &atom->line[kr];

	  if (line->PRD) {
	    Nlamu = 2*atmos.Nrays * line->Nlambda;
	    
	    line->frac = matrix_double(Nlamu, atmos.Nspace);
	    line->id0  = matrix_int(Nlamu, atmos.Nspace);
	    line->id1  = matrix_int(Nlamu, atmos.Nspace);

	    for (la = 0;  la < line->Nlambda;  la++) {
	      for (mu = 0;  mu < atmos.Nrays;  mu++) {
		for (to_obs = 0;  to_obs <= 1;  to_obs++) {
		  sign = (to_obs) ? 1.0 : -1.0;
		  lamu = 2*(atmos.Nrays*la + mu) + to_obs;

		  for (k = 0;  k < atmos.Nspace;  k++) {

		    /* ---  Wavelength in local rest frame -- ------- */
		    
		    lag=line->lambda[la] *
		      (1.0 + spectrum.v_los[mu][k]*sign/CLIGHT);

		    if (lag <= line->lambda[0]) {
		      
		      /* --- out of the lambda table then
			     constant extrapolation -- ------------- */
		      
		      line->frac[lamu][k] = 0.0;
		      line->id0[lamu][k]  = 0;
		      line->id1[lamu][k]  = 1;
		      
		    } else if (lag >= line->lambda[line->Nlambda-1]) {
		      
		      /* --- Out of the lambda table then
			     constant extrapolation -- ------------- */
		      
		      line->frac[lamu][k] = 1.0;
		      line->id0[lamu][k]  = line->Nlambda-2;
		      line->id1[lamu][k]  = line->Nlambda-1;
		      
		    } else {
		      
		      /* --- Locate index of line->lambda of point
			     directly to the left of lag --- ------- */
		      
		      Locate(line->Nlambda,line->lambda,lag,&ilow);
		      line->frac[lamu][k] = (lag-line->lambda[ilow]) /
			(line->lambda[ilow+1] - line->lambda[ilow]);
		      line->id0[lamu][k] = ilow;
		      line->id1[lamu][k] = ilow+1;
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }

    /* --- Precompute Jgas interpolation coefficients if requested -- */
    
    if (!input.prdh_limit_mem) {
      lambda = spectrum.lambda;

      /* --- keeps track of where to get indices and interpolation
             coefficients in spectrum.iprhh and spectrum.cprdh -- --- */
      
      spectrum.nc =
	(int *) malloc(2*atmos.Nrays*spectrum.Nspect*atmos.Nspace *
		       sizeof(int));

      for (la = 0;  la < spectrum.Nspect;  la++) {
	for (mu = 0;  mu < atmos.Nrays;  mu++) {
	  for (to_obs = 0;  to_obs <= 1;  to_obs++) {

	    sign = (to_obs) ? 1.0 : -1.0;

	    for (k = 0;  k < atmos.Nspace;  k++) {

	      lamuk = la * (atmos.Nrays*2*atmos.Nspace)
		+ mu     * (2*atmos.Nspace)
		+ to_obs * (atmos.Nspace)
		+ k ;

	      ncoef = 0;

	      /* --- Previous, current and next wavelength shifted to
		     gas rest frame --                 -------------- */
	      
	      fac = (1.0 + spectrum.v_los[mu][k] * sign / CLIGHT);
	      lambda_prv = lambda[ MAX(la-1,0)                 ] * fac;
	      lambda_gas = lambda[ la                          ] * fac;
	      lambda_nxt = lambda[ MIN(la+1,spectrum.Nspect-1) ] * fac;

	      /* --- Do lambda_prv and lambda_gas bracket
		     lambda points? --                 -------------- */
	      
	      if (lambda_prv !=  lambda_gas) {
		dl= lambda_gas - lambda_prv;
		for (idx = 0; idx < spectrum.Nspect ; idx++) {
		  if (lambda[idx] > lambda_prv &&
		       lambda[idx] <= lambda_gas) ncoef = ncoef+1;
		}
	      } else {
		/* --- edge case, use constant extrapolation for
		       lambda[idx]<lambda gas --       -------------- */
		
		for (idx = 0; idx < spectrum.Nspect ; idx++) {
		  if (lambda[idx] <=  lambda_gas) ncoef=ncoef+1;
		}
	      }
	      /* --- Do lambda_gas and lambda_nxt bracket
		     lambda points? --                 -------------- */
	      
	      if (lambda_gas != lambda_nxt) {
		dl= lambda_nxt - lambda_gas;
		for (idx = 0; idx < spectrum.Nspect ; idx++) {
		  if (lambda[idx] > lambda_gas &&
		      lambda[idx] < lambda_nxt) ncoef=ncoef+1;
		}
	      } else {
		/* --- Edge case, use constant extrapolation for
		       lambda[idx]>lambda gas --       -------------- */
		
		for (idx = 0; idx < spectrum.Nspect ; idx++) {
		  if (lambda[idx] >=  lambda_gas) ncoef = ncoef + 1;
		}
	      }
	      /* --- number of point this lambda contributes to is
	 	     computed as a difference --       -------------- */
	      
	      if (lamuk == 0) {
		spectrum.nc[lamuk] = ncoef;
	      } else {
		spectrum.nc[lamuk] = spectrum.nc[lamuk-1] + ncoef;
	      }
	    }
	  }
	}
      }
      /* --- Now we know the number of interpolation coefficients,
             it's stored in the last element of spectrum.nc,
	     so allocate space --                      -------------- */
      
      idx = spectrum.nc[2*atmos.Nrays * spectrum.Nspect *
			atmos.Nspace - 1];
      spectrum.iprdh = (int *)    malloc( idx * sizeof(int   ));
      spectrum.cprdh = (double *) malloc( idx * sizeof(double));

      /* --- Run through all lamuk points again, and now store indices
             to lambda array and the corresponding interpolation
             coefficients --                           -------------- */
      
      for (la = 0;  la < spectrum.Nspect;  la++) {
	for (mu = 0;  mu < atmos.Nrays;  mu++) {
	  for (to_obs = 0;  to_obs <= 1;  to_obs++) {

	    sign = (to_obs) ? 1.0 : -1.0;

	    for (k = 0;  k < atmos.Nspace;  k++) {

	      lamuk = la * (atmos.Nrays*2*atmos.Nspace)
		+ mu     * (2*atmos.Nspace)
		+ to_obs * (atmos.Nspace)
		+ k ;

	      /* --- Starting index for storage for this lamuk point- */
	      
	      lc = (lamuk==0) ? 0 : spectrum.nc[lamuk-1];

	      /* --- Previous, current and next wavelength shifted
		     to gas rest frame --              -------------- */
	      
	      fac = (1.0 + spectrum.v_los[mu][k] * sign/CLIGHT);
	      lambda_prv = lambda[ MAX(la-1,0)                 ] * fac;
	      lambda_gas = lambda[ la                          ] * fac;
	      lambda_nxt = lambda[ MIN(la+1,spectrum.Nspect-1) ] * fac;

	      /* --- Do lambda_prv and lambda_gas bracket
		     lambda points? --                --------------- */
	      
	      if (lambda_prv !=  lambda_gas) {
		dl = lambda_gas - lambda_prv;
		for (idx = 0;  idx < spectrum.Nspect;  idx++) {
		  if (lambda[idx] > lambda_prv &&
		      lambda[idx] <= lambda_gas) {
		    
		    /* --- Bracketed point found --   --------------- */
		    
		    spectrum.iprdh[lc] = idx;
		    spectrum.cprdh[lc] = (lambda[idx]-lambda_prv) / dl;
		    lc++;
		  }
		}
	      } else {
		/* --- Edge case, use constant extrapolation for
		       lambda[idx]<lambda gas --      --------------- */
		
		for (idx = 0;  idx < spectrum.Nspect;  idx++) {
		  if (lambda[idx] <=  lambda_gas)  {
		    spectrum.iprdh[lc] = idx;
		    spectrum.cprdh[lc] = 1.0;
		    lc++;
		  }
		}
	      }
	      /* --- Do lambda_gas and lambda_nxt bracket
		     lambda points? --                 -------------- */
	      
	      if (lambda_gas != lambda_nxt) {
		dl= lambda_nxt - lambda_gas;
		for (idx = 0;  idx < spectrum.Nspect;  idx++) {
		  if (lambda[idx] > lambda_gas &&
		      lambda[idx] < lambda_nxt) {
		    
		    /* --- Bracketed point found --    -------------- */
		    
		    spectrum.iprdh[lc] = idx;
		    spectrum.cprdh[lc] = 1.0 -
		      (lambda[idx] - lambda_gas) / dl;
		    lc++;
		  }
		}
	      } else {
		/* --- Edge case, use constant extrapolation for
		       lambda[idx]>lambda gas --       -------------- */
		
		for (idx = 0;  idx < spectrum.Nspect;  idx++) {
		  if (lambda[idx] >=  lambda_gas)  {
		    spectrum.iprdh[lc] = idx;
		    spectrum.cprdh[lc] = 1.0;
		    lc++;
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }

  /* --- Allocate space for the emergent intensity --  -------------- */

  switch (topology) {
  case ONE_D_PLANE:
    spectrum.I = matrix_double(spectrum.Nspect, atmos.Nrays);
    if (atmos.Stokes || input.backgr_pol) {
      spectrum.Stokes_Q = matrix_double(spectrum.Nspect, atmos.Nrays);
      spectrum.Stokes_U = matrix_double(spectrum.Nspect, atmos.Nrays);
      spectrum.Stokes_V = matrix_double(spectrum.Nspect, atmos.Nrays);
    }
    break;
  case TWO_D_PLANE:
    Nsr = spectrum.Nspect * atmos.Nrays;
    spectrum.I = matrix_double(Nsr, atmos.N[0]);
    if (atmos.Stokes || input.backgr_pol) {
      spectrum.Stokes_Q = matrix_double(Nsr, atmos.N[0]);
      spectrum.Stokes_U = matrix_double(Nsr, atmos.N[0]);
      spectrum.Stokes_V = matrix_double(Nsr, atmos.N[0]);
    }
    break;
  case THREE_D_PLANE:
    spectrum.I = matrix_double(spectrum.Nspect * atmos.Nrays, 
			       atmos.N[0] * atmos.N[1]);
    if (atmos.Stokes || input.backgr_pol) {
      Nsr    = spectrum.Nspect * atmos.Nrays;
      Nplane = atmos.N[0] * atmos.N[1];

      spectrum.I = matrix_double(Nsr, Nplane);
      if (atmos.Stokes || input.backgr_pol) {
	spectrum.Stokes_Q = matrix_double(Nsr, Nplane);
	spectrum.Stokes_U = matrix_double(Nsr, Nplane);
	spectrum.Stokes_V = matrix_double(Nsr, Nplane);
      }
    }
    break;
  case SPHERICAL_SYMMETRIC:
    spectrum.I = matrix_double(spectrum.Nspect, atmos.Nrays);
    if (atmos.Stokes) {
      Error(ERROR_LEVEL_2, routineName,
	    "Cannot do a full Stokes solution in spherical geometry");
    }    
    break;
  default:
    sprintf(messageStr, "Unknown topology (%d)", topology);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }
  /* --- Read angle-averaged intensity from previous run if necessary,
         and open file for J in case option for limited memory is set */

  spectrum.fd_J   = -1;
  spectrum.fd_J20 = -1;
  oflag = 0;
  openJfile = FALSE;

  if (input.startJ == OLD_J) {
    if (spectrum.updateJ) {
      strcpy(permission, "r+"); 
      oflag |= O_RDWR;
    } else {
      strcpy(permission, "r"); 
      oflag |= O_RDONLY;
    }
    openJfile = TRUE;
  } else {
    if (input.limit_memory) {
      strcpy(permission, "w+"); 
      oflag |= (O_RDWR | O_CREAT);
      openJfile = TRUE;
    }
  }
  if (openJfile) {
    if ((spectrum.fd_J = open(input.JFile, oflag, PERMISSIONS)) == -1) {
      sprintf(messageStr,
	      "Unable to open input file %s with permission %s",
	      input.JFile, permission);
      Error(ERROR_LEVEL_2, routineName, messageStr);
    }
    if (input.backgr_pol) {
      if ((spectrum.fd_J20 = open(J20_DOT_OUT, oflag,
				  PERMISSIONS)) == -1) {
	sprintf(messageStr,
		"Unable to open input file %s with permission %s",
		J20_DOT_OUT, permission);
	Error(ERROR_LEVEL_2, routineName, messageStr);
      }
    }
  }
  if (input.limit_memory) {
    if (oflag & O_CREAT) {
      J = (double *) malloc(atmos.Nspace * sizeof(double));

      /* --- Initialize J file with zeroes --          -------------- */
      
      for (k = 0;  k < atmos.Nspace;  k++) J[k] = 0.0;
      for (nspect = 0;  nspect < spectrum.Nspect;  nspect++)
	writeJlambda(nspect, J);

      free(J);

      if (input.backgr_pol) {
	J20 = (double *) malloc(atmos.Nspace * sizeof(double));
	for (k = 0;  k < atmos.Nspace;  k++) J20[k] = 0.0;
	for (nspect = 0;  nspect < spectrum.Nspect;  nspect++)
	  writeJ20lambda(nspect, J20);

	free(J20);
      }
    }
  } else {
    if (input.startJ == OLD_J) {

      /* --- Fill matrix J with old values from previous run ----- -- */

      for (nspect = 0;  nspect < spectrum.Nspect;  nspect++)
	readJlambda(nspect, spectrum.J[nspect]);
    
      close(spectrum.fd_J);
      spectrum.fd_J = -1;

      if (input.backgr_pol) {
	for (nspect = 0;  nspect < spectrum.Nspect;  nspect++)
	  readJ20lambda(nspect, spectrum.J20[nspect]);
    
	close(spectrum.fd_J20);
	spectrum.fd_J20 = -1;
      }
    }
  }

  /* --- Look for Jgas and read, otherwise use spectrum.J --- ------- */
  
  if (atmos.NPRDactive > 0 &&
      input.PRD_angle_dep == PRD_ANGLE_APPROX) {
      
    fp = fopen("Jgas.dat", "r");
    if (fp) {
	
      /* --- File exists --                            -------------- */
	
      fclose(fp);
      readJgas(spectrum.Jgas);
      sprintf(messageStr, "Read spectrum.Jgas from file.");
      Error(MESSAGE, routineName, messageStr);
    } else {
	
      /* --- File does not exist --                    -------------- */
	
      sprintf(messageStr,
	      "Jgas.dat does not exist,setting spectrum.Jgas spectrum.J.");
      Error(WARNING, routineName, messageStr);
	
      for (k = 0;  k < atmos.Nspace;  k++) {
	for (nspect = 0;  nspect < spectrum.Nspect;  nspect++) {
	  spectrum.Jgas[nspect][k] = spectrum.J[nspect][k];
	}
      }
    }
  }

  /* --- Need storage for angle-dependent specific intensities for
         angle-dependent PRD --                        -------------- */

  if (atmos.NPRDactive > 0  &&
      input.PRD_angle_dep == PRD_ANGLE_DEP) {
    
    oflag = 0;
    if (input.startJ == OLD_J) {
      if (spectrum.updateJ) {
	strcpy(permission, "r+");
	oflag |= O_RDWR;
      } else {
	strcpy(permission, "r");
	oflag |= O_RDONLY;
      }
    } else {
      strcpy(permission, "w+");
      oflag |= (O_RDWR | O_CREAT);
    }
    if ((spectrum.fd_Imu = open(IMU_FILENAME, oflag, PERMISSIONS)) == -1) {
      sprintf(messageStr, "Unable to open %s file %s with permission %s",
	      (spectrum.updateJ) ? "update" : "input",
	      IMU_FILENAME, permission);
      Error(ERROR_LEVEL_2, routineName, messageStr);
    }
    /* --- Fill the index list that keeps track of the location
           of intensity Imu in file spectrum.fd_Imu at wavelength
           corresponding to nspect. --                 -------------- */

    spectrum.PRDindex = (int *) malloc(spectrum.Nspect * sizeof(int));
    index = 0;
    for (nspect = 0;  nspect < spectrum.Nspect;  nspect++) {
      if (containsPRDline(&spectrum.as[nspect])) {
	spectrum.PRDindex[nspect] = index;
        index++;
      }
    }
  }
  for (nact = 0;  nact < atmos.Nactiveatom;  nact++) {
    atom = atmos.activeatoms[nact];

    /* --- Allocate memory for the rate equation matrix -- ---------- */

    atom->Gamma = matrix_double(SQ(atom->Nlevel), atmos.Nspace);

    /* --- Initialize the mutex lock for the operator Gamma if there
           are more than one threads --                -------------- */

    if (input.Nthreads > 0) {
      if ((status = pthread_mutex_init(&atom->Gamma_lock, NULL))) {
	sprintf(messageStr, "Unable to initialize mutex_lock, status = %d",
		status);
	Error(ERROR_LEVEL_2, routineName, messageStr);
      }
    }

    switch(atom->initial_solution) {
    case LTE_POPULATIONS:
      for (i = 0;  i < atom->Nlevel;  i++) {
	for (k = 0;  k < atmos.Nspace;  k++)
	  atom->n[i][k] = atom->nstar[i][k];
      }
      break;

    case ZERO_RADIATION:
      hc_k   = (HPLANCK * CLIGHT) / (KBOLTZMANN * NM_TO_M);
      twoc   = 2.0*CLIGHT / CUBE(NM_TO_M);
      fourPI = 4.0 * PI;

      initGammaAtom(atom, 1);

      /* --- Then add radiative contributions of active transitions --  */

      for (nspect = 0;  nspect < spectrum.Nspect;  nspect++) {
	as = spectrum.as + nspect;

	for (n = 0;  n < as->Nactiveatomrt[nact];  n++) {
	  switch (as->art[nact][n].type) {
	  case ATOMIC_LINE:
	    line = as->art[nact][n].ptype.line;
	    la = nspect - line->Nblue;
	    i  = line->i;
	    j  = line->j;
	    ij = i*atom->Nlevel + j;
	    
	    if (la == 0) {
	      for (k = 0;  k < atmos.Nspace;  k++)
		atom->Gamma[ij][k] += line->Aji;
	    }
	    break;

	  case ATOMIC_CONTINUUM:
	    continuum = as->art[nact][n].ptype.continuum;
	    la = nspect - continuum->Nblue;
	    i  = continuum->i;
	    j  = continuum->j;
	    ij = i*atom->Nlevel + j;

	    wla = fourPI * getwlambda_cont(continuum, la) /
	      continuum->lambda[la];
	    twohnu3_c2 = twoc / CUBE(continuum->lambda[la]);
	    for (k = 0;  k < atmos.Nspace;  k++) {
	      gijk = atom->nstar[i][k]/atom->nstar[j][k] *
		exp(-hc_k/(continuum->lambda[la] * atmos.T[k]));
	      atom->Gamma[ij][k] += gijk * twohnu3_c2 *
		continuum->alpha[la]*wla;
	    }
	    break;
	  default:
	    break;
	  }
	}
      }
      /* --- Solve statistical equilibrium equations --  ------------ */

      statEquil(atom, (input.isum == -1) ? 0 : input.isum);
      break;

    case OLD_POPULATIONS:
      readPopulations(atom);
      break;

    default:;
    break;
    }
  }
  /* --- Now the molecules that are active --          -------------- */
  
  for (nact = 0;  nact < atmos.Nactivemol;  nact++) {
    molecule = atmos.activemols[nact];

    /* --- Calculate the LTE vibration level populations here. They
           cannot be calculated yet in readMolecule since chemical
           equilibrium has to be established first --  -------------- */

    for (i = 0;  i < molecule->Nv;  i++) {
      for (k = 0;  k < atmos.Nspace;  k++)
	molecule->nvstar[i][k] = molecule->n[k] *
	  molecule->pfv[i][k] / molecule->pf[k];
    }
    /* --- Allocate memory for the rate equation matrix -- ---------- */

    molecule->Gamma = matrix_double(SQ(molecule->Nv), atmos.Nspace);

    /* --- Initialize the mutex lock for the operator Gamma if there
           are more than one thread --                 -------------- */

    if (input.Nthreads > 0) {
      if ((status = pthread_mutex_init(&molecule->Gamma_lock, NULL))) {
	sprintf(messageStr, "Unable to initialize mutex_lock, status = %d",
		status);
	Error(ERROR_LEVEL_2, routineName, messageStr);
      }
    }

    switch(molecule->initial_solution) {

    case LTE_POPULATIONS:
      for (i = 0;  i < molecule->Nv;  i++) {
	for (k = 0;  k < atmos.Nspace;  k++)
	  molecule->nv[i][k] = molecule->nvstar[i][k];
      }
      break;
      
    case OLD_POPULATIONS:
      readMolPops(molecule);
      break;

    default:;
    }

    /* --- Calculate collisions for molecule (must be done here because
           rotation-vibration transitions are dominated by hydrogen and
           H2 collisions for which chemical equilibrium needs to be
           established first --                        -------------- */

    if (strstr(molecule->ID, "CO"))
      COcollisions(molecule);
    else if (strstr(molecule->ID, "H2"))
      H2collisions(molecule);
    else {
      sprintf(messageStr, "Collisions for molecule %s not implemented\n",
	      molecule->ID);
      Error(ERROR_LEVEL_2, routineName, messageStr);
    }
  }
}
/* ------- end ---------------------------- initSolution.c ---------- */
