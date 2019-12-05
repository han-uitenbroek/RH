#include <stdlib.h>
#include <math.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "spectrum.h"
#include "geometry.h"
#include "accelerate.h"
#include "inputs.h"
#include "error.h"
#include "statistics.h"
#include "constant.h"
#include "xdr.h"

enum  Topology topology = TWO_D_PLANE;

Atmosphere atmos;
Geometry geometry;
Spectrum spectrum;
ProgramStats stats;
CommandLine commandline;
InputData input;

char messageStr[MAX_LINE_SIZE];

struct Ng *Ng_S;

   int  Nlambda;
double *Bp, *chi, *S, **Iemerge, *epsilon, *phi, *wlamb, wphi;
  FILE *fpin, *fpout;
static  XDR  xdrs;

#define N_INPUTS  8

int main(int argc, char *argv[])
{
  register int k, l, n, la;

  int     Nspace, result, NmaxIter, Ngdelay, Ngorder, Ngperiod, btype[3],
          inputs[N_INPUTS], nwrite, Nx, Nz;
  double  iterLimit, *lambda, Adamp;

  stats.printCPU = TRUE;
  commandline.quiet = FALSE;
  commandline.logfile = stderr;
  input.Eddington = FALSE;

  input.S_interpolation = S_PARABOLIC;

  getCPU(0, TIME_START, NULL);
  SetFPEtraps();

  fpin  = fopen("2dinput.dat", "r");
  xdrstdio_create(&xdrs, fpin, XDR_DECODE);

  result = xdr_vector(&xdrs, (char *) inputs, N_INPUTS,
		      sizeof(int), (xdrproc_t) xdr_int);

  atmos.angleSet.set = (enum angleset) inputs[0];
  if (atmos.angleSet.set == SET_GL) {
     result = xdr_int(&xdrs, &atmos.angleSet.Ninclination);
     result = xdr_int(&xdrs, &atmos.angleSet.Nazimuth);
  }
  result = xdr_double(&xdrs, &iterLimit);
  result = xdr_double(&xdrs, &Adamp);

  Nx = geometry.Nx = inputs[1];
  Nz = geometry.Nz = inputs[2];
  NmaxIter = inputs[3];
  Ngdelay  = inputs[4];  Ngorder = inputs[5];  Ngperiod = inputs[6];
  Nlambda  = inputs[7];
  Nspace   = atmos.Nspace = geometry.Nx * geometry.Nz;

  result = xdr_vector(&xdrs, (char *) btype, 3,
		      sizeof(int), (xdrproc_t) xdr_int);
  geometry.hboundary = (enum boundary) btype[0];
  for (n = 0;  n < 2;  n++)
    geometry.bvalue[n] = (enum boundval) btype[1+n];

  /* --- Check the validity of boundary conditions and values -- ---- */ 

  switch (geometry.hboundary) {
  case FIXED:     break;
  case PERIODIC:  break;
  default:
    sprintf(messageStr, "Invalid horizontal boundary condition: %d",
	    geometry.hboundary);
    Error(ERROR_LEVEL_2, argv[0], messageStr);
    break;
  }

  switch (geometry.bvalue[TOP]) {
  case IRRADIATED:   break;
  case ZERO:         break;
  case THERMALIZED:  break;
  default:
    sprintf(messageStr, "Invalid boundary value at TOP: %d\n",
	    geometry.bvalue[TOP]);
    Error(ERROR_LEVEL_2, argv[0], messageStr);
    break;
  }
  switch (geometry.bvalue[BOTTOM]) {
  case IRRADIATED:   break;
  case ZERO:         break;
  case THERMALIZED:  break;
  default:
    sprintf(messageStr, "Invalid boundary value at BOTTOM: %d",
	    geometry.bvalue[BOTTOM]);
    Error(ERROR_LEVEL_2, argv[0], messageStr);
    break;
  }

  /* --- Get increments in x, store and check for monotonicity -- --- */

  geometry.dx = (double *) malloc(Nx * sizeof(double));
  result = xdr_vector(&xdrs, (char *) geometry.dx, Nx,
		       sizeof(double), (xdrproc_t) xdr_double);
  for (l = 0;  l < ((geometry.hboundary == PERIODIC) ? Nx : Nx-1);  l++) {
    geometry.dx[l] *= KM_TO_M;
    if (geometry.dx[l] <= 0.0) {
      sprintf(messageStr, "At l = %d:\n x does not increase strictly "
              "monotonically towards the right", l);
      Error(ERROR_LEVEL_2, argv[0], messageStr);
    } 
  }
  geometry.x = (double *) malloc(Nx * sizeof(double));
  geometry.x[0] = 0.0;
  for (l = 0;  l < Nx-1;  l++)
    geometry.x[l+1] = geometry.x[l] + geometry.dx[l];

  /* --- Get vertical grid --                          -------------- */

  geometry.z = (double *) malloc(Nz * sizeof(double));
  result = xdr_vector(&xdrs, (char *) geometry.z, Nz,
		      sizeof(double), (xdrproc_t) xdr_double);
  for (k = 0;  k < Nz;  k++) geometry.z[k] *= KM_TO_M;

  geometry.dz = (double *) malloc(Nz * sizeof(double));
  for (k = 0;  k < Nz-1;  k++) {
    geometry.dz[k] = geometry.z[k] - geometry.z[k+1];
    if (geometry.dz[k] <= 0.0) {
      sprintf(messageStr, "At k = %d:\n z does not decrease strictly "
              "monotonically towards the bottom", k);
      Error(ERROR_LEVEL_2, argv[0], messageStr);
    }
  }
  geometry.dz[Nz-1] = 0.0;


  chi = (double *) malloc(Nspace * sizeof(double));
  S   = (double *) malloc(Nspace * sizeof(double));
  Bp  = (double *) malloc(Nspace * sizeof(double));
  epsilon = (double *) malloc(Nspace * sizeof(double));

  result = xdr_vector(&xdrs, (char *) chi, Nspace,
		       sizeof(double), (xdrproc_t) xdr_double);
  result = xdr_vector(&xdrs, (char *) Bp, Nspace,
		       sizeof(double), (xdrproc_t) xdr_double);
  result = xdr_vector(&xdrs, (char *) epsilon, Nspace,
		       sizeof(double), (xdrproc_t) xdr_double);

  getBoundary(&atmos, &geometry);
  getAngleQuadr(&geometry);
  atmos.Nrays = geometry.Nrays;
  atmos.wmu   = geometry.wmu;
  fillMesh(&geometry);

  lambda = (double *) malloc(Nlambda * sizeof(double));
  phi    = (double *) malloc(Nlambda * sizeof(double));
  wlamb  = (double *) malloc(Nlambda * sizeof(double));

  result = xdr_vector(&xdrs, (char *) lambda, Nlambda,
		      sizeof(double), (xdrproc_t) xdr_double);
  wlamb[0] = (lambda[1] - lambda[0]);
  for (la = 1;  la < Nlambda-1;  la++)
    wlamb[la] = (lambda[la+1] - lambda[la-1]);
  wlamb[Nlambda-1] = (lambda[Nlambda-1] - lambda[Nlambda-2]);

  wphi = 0.0;
  for (la = 0;  la < Nlambda;  la++) {
    phi[la] = Voigt(Adamp, lambda[la], NULL, RYBICKI)/SQRTPI;
    wphi   += wlamb[la]*phi[la];
  }
  wphi = 1.0/wphi;

  xdr_destroy(&xdrs);
  fclose(fpin);

  for (k = 0;  k < Nspace;  k++)  S[k] = Bp[k];
  Ng_S = NgInit(Nspace, Ngdelay, Ngorder, Ngperiod, S);

  Iterate(NmaxIter, iterLimit);

  nwrite = fwrite(S, sizeof(double), Nspace, stdout);
  printTotalCPU();
}
/* ------- end ---------------------------- solve2d.c --------------- */

/* ------- begin -------------------------- getBoundary.c ----------- */

void getBoundary(Atmosphere *atmos, Geometry *geometry)
{
  int  result, Nspect = Nlambda;

  /* --- Read boundary conditions and irradiation. First VERTICAL
         boundary conditions --                            --------- */

  switch (geometry->bvalue[TOP]) {
  case IRRADIATED:
    geometry->Itop = matrix_double(Nspect, geometry->Nx);
    result = xdr_vector(&xdrs, (char *) geometry->Itop[0],
			Nspect*geometry->Nx,
			sizeof(double), (xdrproc_t) xdr_double);
    break;
  default:
    geometry->Itop = NULL;
    break;
  }

  switch (geometry->bvalue[BOTTOM]) {
  case IRRADIATED:
    geometry->Ibottom = matrix_double(Nspect, geometry->Nx);
    result = xdr_vector(&xdrs, (char *) geometry->Ibottom[0],
			Nspect*geometry->Nx,
			sizeof(double), (xdrproc_t) xdr_double);
    break;
  default:
    geometry->Ibottom = NULL;
    break;
  }

  /* --- Get the HORIZONTAL boundary conditions --        ----------- */

  switch (geometry->hboundary) {
  case FIXED:
    geometry->Ileft = matrix_double(Nspect, geometry->Nz);
    result = xdr_vector(&xdrs, (char *) geometry->Ileft[0],
			Nspect*geometry->Nz,
                        sizeof(double), (xdrproc_t) xdr_double);
    geometry->Iright = matrix_double(Nspect, geometry->Nz);
    result = xdr_vector(&xdrs, (char *) geometry->Iright[0],
			Nspect*geometry->Nz,
                        sizeof(double), (xdrproc_t) xdr_double);
    break;
  case PERIODIC:
    geometry->Ileft  = NULL;
    geometry->Iright = NULL;
    break;
  default:
    break;
  }
}
/* ------- end ---------------------------- getBoundary.c ----------- */
