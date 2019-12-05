/* ------- file: -------------------------- chartflow.c -------------

       Version:       rh2.0, 2-D Cartesian
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Fri May 25 14:00:54 2018 --

       --------------------------                      ----------RH-- */

/* --- Compute the x- and z components of the flux vector on the
       whole grid for given wavelength indices nspect1  ..... nspectN.

       The equation of transfer is solved with the shortcharacteristics
       method.

  See: P. B. Kunasz and L. H. Auer 1988, JQSRT 39, p. 67


       Expects input file ``flow.input'' containing one line of the form

         Nspect  wave_index1  ....   wave_indexNspect
       --                                              -------------- */

#include <fcntl.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "geometry.h"
#include "spectrum.h"
#include "background.h"
#include "statistics.h"
#include "inputs.h"
#include "error.h"
#include "xdr.h"

#define COMMENT_CHAR    "#"
#define FLOW_DOT_INPUT  "flow.input"


/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

enum Topology topology = TWO_D_PLANE;

Atmosphere atmos;
Geometry geometry;
Spectrum spectrum;
ProgramStats stats;
InputData input;
CommandLine commandline;
char messageStr[MAX_LINE_SIZE];


/* ------- begin -------------------------- chartflow.c ------------- */

int main(int argc, char *argv[])
{
  register int k, mu, n, m;

  char    chartFileName[10], inputLine[MAX_LINE_SIZE];
  bool_t  result, exit_on_EOF, crosscoupling, to_obs, initialize,
          boundbound, solve_Stokes, PRD_angle_dep, polarized_as,
          polarized_c, angle_dep, analyze_output, equilibria_only;
  int     Nspect, Nread, *wave_index, checkPoint;
  double  wmux, wmuz, *I, *J, *Psi, *S, *chi, *chi_as, *eta_as,
         *Fx, *Fz;
  FILE   *fp_out = NULL, *fp_flow = NULL;
  XDR     xdrs;
  Atom *atom;
  Molecule *molecule;
  ActiveSet *as;

  setOptions(argc, argv);
  getCPU(0, TIME_START, NULL);
  SetFPEtraps();

  /* --- Read input data and initialize --             -------------- */

  readInput();
  spectrum.updateJ = FALSE;

  getCPU(1, TIME_START, NULL);
  readAtmos(&atmos, &geometry);
  if (atmos.Stokes) Bproject();
  fillMesh(&geometry);
  
  input.startJ = OLD_J;

  readAtomicModels();
  readMolecularModels();
  SortLambda();

  getBoundary(&atmos, &geometry);

  /* --- Solve chemical equilibrium and open background opacity file- */
  
  Background(analyze_output=FALSE, equilibria_only=TRUE);

  if ((atmos.fd_background =
       open(input.background_File, O_RDONLY, 0)) == -1) {
    sprintf(messageStr, "Unable to open inputfile %s",
	    input.background_File);
    Error(ERROR_LEVEL_2, argv[0], messageStr);
  }
  readBRS();

  getProfiles();
  initSolution();
  initScatter();

  /* --- Get wavelength indices for flow field --      -------------- */

  if ((fp_flow = fopen(FLOW_DOT_INPUT, "r")) == NULL) {
    sprintf(messageStr, "Unable to open inputfile %s", FLOW_DOT_INPUT);
    Error(ERROR_LEVEL_2, argv[0], messageStr);
  }

  getLine(fp_flow, COMMENT_CHAR, inputLine, exit_on_EOF=TRUE);
  Nread = sscanf(strtok(inputLine, " "), "%d", &Nspect);
  wave_index = (int *) malloc(Nspect * sizeof(int));
  for (n = 0;  n < Nspect;  n++)
    Nread += sscanf(strtok(NULL, " "), "%d", wave_index+n);
  checkNread(Nread, Nspect+1, argv[0], checkPoint=2);
  fclose(fp_flow);
  

  I      = (double *) malloc(atmos.Nspace * sizeof(double));
  Psi    = (double *) malloc(atmos.Nspace * sizeof(double));
  chi    = (double *) malloc(atmos.Nspace * sizeof(double));
  chi_as = (double *) malloc(atmos.Nspace * sizeof(double));
  eta_as = (double *) malloc(atmos.Nspace * sizeof(double));
  S      = (double *) malloc(atmos.Nspace * sizeof(double));
  
  Fx = (double *) malloc(atmos.Nspace * sizeof(double));
  Fz = (double *) malloc(atmos.Nspace * sizeof(double));

  /* --- Solve radiative transfer equation at specified wavelengths - */
  
  getCPU(1, TIME_POLL, "Total initialize");

  for (n = 0;  n < Nspect;  n++) {
    if (wave_index[n] < 0  ||  wave_index[n] >= spectrum.Nspect) {
      sprintf(messageStr, "Illegal value of wave_index[n]: %4d\n"
	      "Value has to be between 0 and %4d\n", 
	      wave_index[n], spectrum.Nspect);
      Error(ERROR_LEVEL_2, argv[0], messageStr);
      continue;
    }
    sprintf(messageStr, "Processing n = %4d, lambda = %9.3f [nm]\n",
	    wave_index[n], spectrum.lambda[wave_index[n]]);
    Error(MESSAGE, NULL, messageStr);

    /* --- Retrieve active set of transitions as and read current
           mean intensities at wavelength nspect. --   -------------- */

    as = &spectrum.as[wave_index[n]];
    alloc_as(wave_index[n], crosscoupling=FALSE);

    if (input.limit_memory) {
      J = (double *) malloc(atmos.Nspace * sizeof(double));
      readJlambda(wave_index[n], J);
    } else {
      J = spectrum.J[wave_index[n]];
    }

    /* --- Check whether current active set includes a bound-bound
           and/or polarized transition and/or angledependent PRD
           transition. Otherwise, only angle-independent opacity and
           source functions are needed --              -------------- */ 


    boundbound    = containsBoundBound(as);
    PRD_angle_dep = (containsPRDline(as) && input.PRD_angle_dep);
    polarized_as  = containsPolarized(as);
    polarized_c   = atmos.backgrflags[wave_index[n]].ispolarized;
    solve_Stokes  = (input.StokesMode == FULL_STOKES &&
		     (polarized_as || polarized_c));
    angle_dep     = (polarized_as || polarized_c || PRD_angle_dep ||
		     (atmos.moving &&
	     (boundbound || atmos.backgrflags[wave_index[n]].hasline)));


    if (angle_dep) {
      for (mu = 0;  mu < geometry.Nrays;  mu++) {
	wmux = geometry.mux[mu] * geometry.wmu[mu];
	wmuz = geometry.muz[mu] * geometry.wmu[mu];
	for (to_obs = 0;  to_obs <= 1;  to_obs++) {
	  initialize = (mu == 0 && to_obs == 0);

	  if (initialize || atmos.backgrflags[wave_index[n]].hasline)
	    readBackground(wave_index[n], mu, to_obs);
	  
	  if (initialize || boundbound)
	    Opacity(wave_index[n], mu, to_obs, initialize);

	  for (k = 0;  k < atmos.Nspace;  k++) {
	    chi[k] = as->chi[k] + as->chi_c[k];
	    S[k]   = (as->eta[k] + as->eta_c[k] +
		      as->sca_c[k]*J[k]) / chi[k];
	  }
	  Piecewise_2D(&geometry, wave_index[n], mu, to_obs,
		       chi, S, I, Psi);

	  if (to_obs) {
	    for (k = 0;  k < atmos.Nspace;  k++) {
	      Fx[k] += I[k] * wmux;
	      Fz[k] += I[k] * wmuz;
	    }
	  } else {
	    for (k = 0;  k < atmos.Nspace;  k++) {
	      Fx[k] -= I[k] * wmux;
	      Fz[k] -= I[k] * wmuz;
	    }
	  }
	}
      }
    } else {
      Opacity(wave_index[n], 0, 0, initialize=TRUE);
      readBackground(wave_index[n], 0, 0);
      
      for (k = 0;  k < atmos.Nspace;  k++) {
	chi[k] = chi_as[k] + as->chi_c[k];
	S[k]   = (eta_as[k] + as->eta_c[k] +
		  as->sca_c[k]*J[k]) / chi[k];
      }
      for (mu = 0;  mu < geometry.Nrays;  mu++) {
	wmux = geometry.mux[mu] * geometry.wmu[mu];
	wmuz = geometry.muz[mu] * geometry.wmu[mu];
	for (to_obs = 0;  to_obs <= 1;  to_obs++) {
	  Piecewise_2D(&geometry, wave_index[n], mu,
		       to_obs, chi, S, I, Psi);

	  if (to_obs) {
	    for (k = 0;  k < atmos.Nspace;  k++) {
	      Fx[k] += I[k] * wmux;
	      Fz[k] += I[k] * wmuz;
	    }
	  } else {
	    for (k = 0;  k < atmos.Nspace;  k++) {
	      Fx[k] -= I[k] * wmux;
	      Fz[k] -= I[k] * wmuz;
	    }
	  }
	}
      }
    }
    free_as(wave_index[n], crosscoupling=FALSE);

    /* --- Write the radiation flow to output file --  -------------- */
    
    sprintf(chartFileName, "flow_%4.4d", wave_index[n]);
    if ((fp_out = fopen(chartFileName, "w" )) == NULL) {
      sprintf(messageStr, "Unable to open output file %s", chartFileName);
      Error(ERROR_LEVEL_2, argv[0], messageStr);
    }
    xdrstdio_create(&xdrs, fp_out, XDR_ENCODE);
    
    result = xdr_int(&xdrs, &wave_index[n]);
    result = xdr_vector(&xdrs, (char *) Fx, atmos.Nspace,
			sizeof(double), (xdrproc_t) xdr_double);
    result = xdr_vector(&xdrs, (char *) Fz, atmos.Nspace,
			sizeof(double), (xdrproc_t) xdr_double);

    if (input.limit_memory) free(J);
  }
  xdr_destroy(&xdrs);
  fclose(fp_out);
  printTotalCPU();
}
/* ------- end ---------------------------- chartflow.c ------------- */
