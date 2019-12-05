/* ------- file: -------------------------- formal2d.c -------------- */

/* Formal solution in 2-D geometry of just one ray.
 *
 * Han Uitenbroek
 * Last modified: Fri May 25 13:59:45 2018 --
 */
 
#include <stdlib.h>
#include <math.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "spectrum.h"
#include "geometry.h"
#include "inputs.h"
#include "error.h"
#include "statistics.h"
#include "constant.h"

/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

enum Topology topology = TWO_D_PLANE;

Atmosphere atmos;
Spectrum spectrum;
Geometry geometry;
ProgramStats stats;
CommandLine commandline;
InputData input;
char messageStr[MAX_LINE_SIZE];


/* ------- begin --------------------------             ------------- */

int main(int argc, char *argv[])
{
  register int l, k;

  bool_t  to_observer;
  int     mu, Nread, Nwrite, Nspace, Nspect = 1, nspect, Nx, Nz;
  double *chi, *S, *I, *Psi;

  stats.printCPU = TRUE;
  commandline.quiet = FALSE;
  commandline.logfile = stderr;
  input.Eddington = FALSE;

  input.S_interpolation = S_PARABOLIC;

  getCPU(0, TIME_START, NULL);
  SetFPEtraps();

  Nread = fread(&atmos.angleSet.set, sizeof(int), 1, stdin);
  Nread = fread(&mu, sizeof(int), 1, stdin);
  if (atmos.angleSet.set == SET_GL) {
    Nread = fread(&atmos.angleSet.Ninclination, sizeof(int), 1, stdin);
    Nread = fread(&atmos.angleSet.Nazimuth, sizeof(int), 1, stdin);
  }

  geometry.bvalue[TOP] = geometry.bvalue[BOTTOM] = IRRADIATED;
  Nread = fread(&geometry.hboundary, sizeof(int), 1, stdin);

  Nread = fread(&geometry.Nx, sizeof(int), 1, stdin);
  Nx = geometry.Nx;
  Nread = fread(&geometry.Nz, sizeof(int), 1, stdin);
  Nz = geometry.Nz;
  Nspace = atmos.Nspace = geometry.Nx * geometry.Nz;

  /* --- Get increments in x, store and check for monotonicity -- --- */

  geometry.dx = (double *) malloc(Nx * sizeof(double));
  Nread = fread(geometry.dx, sizeof(double), Nx, stdin);

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
  Nread = fread(geometry.z, sizeof(double), Nz, stdin);
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

  geometry.Itop    = matrix_double(Nspect, geometry.Nx);
  geometry.Ibottom = matrix_double(Nspect, geometry.Nx);
  Nread = fread(geometry.Itop[0], sizeof(double), geometry.Nx, stdin);
  Nread = fread(geometry.Ibottom[0], sizeof(double), geometry.Nx, stdin);

  switch (geometry.hboundary) {
  case FIXED:
    geometry.Ileft = matrix_double(Nspect, geometry.Nz);
    Nread = fread(geometry.Ileft[0], sizeof(double), geometry.Nz, stdin);
    geometry.Iright = matrix_double(Nspect, geometry.Nz);
    Nread = fread(geometry.Iright[0], sizeof(double), geometry.Nz, stdin);
    break;
  case PERIODIC:
    geometry.Ileft  = NULL;
    geometry.Iright = NULL;
    break;
  }

  chi   = (double *) malloc(Nspace * sizeof(double));
  S     = (double *) malloc(Nspace * sizeof(double));
  Nread = fread(chi, sizeof(double), Nspace, stdin);
  Nread = fread(S, sizeof(double), Nspace, stdin);

  I      = (double *) malloc(Nspace * sizeof(double));
  Psi    = (double *) malloc(Nspace * sizeof(double));

  getAngleQuadr(&geometry);
  fillMesh(&geometry);
  getCPU(1, TIME_START, NULL);

  Piecewise_2D(&geometry, nspect=0, mu, to_observer=FALSE, chi, S, I, Psi);

  getCPU(1, TIME_POLL, "Formal Solution down");
  Nwrite = fwrite(I, sizeof(double), Nspace, stdout);
  Nwrite = fwrite(Psi, sizeof(double), Nspace, stdout);

  getCPU(1, TIME_START, NULL);

  Piecewise_2D(&geometry, nspect=0, mu, to_observer=TRUE, chi, S, I, Psi);

  getCPU(1, TIME_POLL, "Formal Solution up");
  Nwrite = fwrite(I, sizeof(double), Nspace, stdout);
  Nwrite = fwrite(Psi, sizeof(double), Nspace, stdout);

  printTotalCPU();
}
/* ------- end ---------------------------- formal2d.c -------------- */
