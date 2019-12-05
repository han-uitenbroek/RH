/* ------- file: -------------------------- initial_xdr.c -----------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Mon Jan 16 17:29:31 2012 --

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

void initSolution(Atom *atom, Molecule *molecule)
{
  const char routineName[] = "initSolution";
  register int k, i, ij, nspect, n, kr, nact;

  char    permission[3];
  bool_t  result, openJfile;
  int     la, j, niter, Nsr, Nplane, index, status, oflag;
  double  gijk, wla, twohnu3_c2, hc_k, twoc, fourPI, *J, *J20;
  ActiveSet *as;
  AtomicLine *line;
  AtomicContinuum *continuum;
  XDR xdrs;

  getCPU(2, TIME_START, NULL);

  /* --- Allocate space for angle-averaged mean intensity -- -------- */

  if (!input.limit_memory) {
    spectrum.J = matrix_double(spectrum.Nspect, atmos.Nspace);

    /* --- If we do background polarization we need space for the
           anisotropy --                               -------------- */

    if (input.backgr_pol)
      spectrum.J20 = matrix_double(spectrum.Nspect, atmos.Nspace);
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
  /* --- Need storage for angle-dependent specific intensities for
         angle-dependent PRD --                        -------------- */

  if (atmos.NPRDactive > 0 && input.PRD_angle_dep) {
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

      initGammaAtom(atom);

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
