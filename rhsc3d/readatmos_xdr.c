/* ------- file: -------------------------- readatmos_xdr.c ---------

       Version:       rh2.0, 3-D Cartesian, short characteristics
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Thu Jun 28 14:42:46 2018 --

       --------------------------                      ----------RH-- */

/* --- Reads the atmospheric data and fills Atmos structure.
       XDR (external data representation) version. --  -------------- */


/* --- Required input:

         Nx, Ny, Nz  -- dimensions of grid             -- 3 4-byte integer
         NHydr       -- number of Hydrogen levels      -- 4-byte integer
                        (this includes the proton number)
         bvalue      -- vertical boundary values       -- 2 4-byte integers
                        IRRADIATED=0, ZERO=1, THERMALIZED=2
                        [TOPvalue, BOTTOMvalue]
         dx          -- horizontal grid dist in x [km] -- double
         dy          -- horizontal grid dist in y [km] -- double
         z           -- vertical scale [km]            -- double[Nz]
         T           -- temperature [K]                -- double[Nx*Ny*Nz]
         ne          -- electron density [m^-3]        -- double[Nx*Ny*Nz]
         vturb       -- microturbulence [km/s]         -- double[Nx*Ny*Nz]
         vx          -- velovity field in x [km/s]     -- double[Nx*Ny*Nz]
         vy          -- velovity field in y [km/s]     -- double[Nx*Ny*Nz]
         vz          -- velovity field in z [km/s]     -- double[Nx*Ny*Nz]
         nH          -- Hydrogen populations [m^-3]    -- double[Nx*Ny*Nz*NHydr]
                        (see background.c in case of LTE hydrogen populations)

       Conventions:

         -- Indexing starts at the TOP of the atmosphere.
         -- z increases from the BOTTOM to the TOP.
         -- Positive vz implies upward motion.

       --                                              -------------- */
 
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "geometry.h"
#include "spectrum.h"
#include "background.h"
#include "constant.h"
#include "statistics.h"
#include "error.h"
#include "inputs.h"
#include "xdr.h"

/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

extern InputData  input;
extern Spectrum   spectrum;
extern char messageStr[];

static XDR xdrs;

/* ------- begin -------------------------- readAtmos.c ------------- */

void readAtmos(Atmosphere *atmos, Geometry *geometry)
{
  const char routineName[] = "readAtmos";
  register int  n, k;

  char   *filename;
  int     result, Nx, Ny, Nplane, Nz, Nspace, bc[2];
  struct  stat statBuffer;

  getCPU(2, TIME_START, NULL);

  /* --- Get abundances of background elements --        ------------ */

  readAbundance(atmos);

  /* --- Open input file for atmospheric model --        ------------ */

  if ((atmos->fp_atmos = fopen(input.atmos_input, "r")) == NULL) {
    sprintf(messageStr, "Unable to open inputfile %s",input.atmos_input);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }
  /* --- Construct atmosID from filename and last modification date - */

  stat(input.atmos_input, &statBuffer);
  if ((filename = strrchr(input.atmos_input, '/')))
    filename++;
  else
    filename = input.atmos_input;
  sprintf(atmos->ID, "%s (%.24s)", filename,
	  asctime(localtime(&statBuffer.st_mtime)));

  xdrstdio_create(&xdrs, atmos->fp_atmos, XDR_DECODE);

  /* --- Read size of horizontal and vertical grids --   ------------ */

  result = xdr_int(&xdrs, &Nx);
  result = xdr_int(&xdrs, &Ny);
  result = xdr_int(&xdrs, &Nz);

  result = xdr_int(&xdrs, &atmos->NHydr);
  if (atmos->NHydr < 2  &&  !atmos->H_LTE) {
    sprintf(messageStr, "NHydr has to be at least 2, not %d to run with"
	    " NLTE hydrogen populations in background\n", atmos->NHydr);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }

  /* --- Keep duplicates of some of the geometrical quantities in
         Atmos structure --                            -------------- */

  atmos->Ndim = 3;
  atmos->N = (int *) malloc(atmos->Ndim * sizeof(int));
  atmos->N[0] = geometry->Nx = Nx;
  atmos->N[1] = geometry->Ny = Ny;
  atmos->N[2] = geometry->Nz = Nz;
  geometry->Nplane = Nplane = Nx * Ny;
  atmos->Nspace    = Nspace = Nplane * Nz;

  /* --- Get the boundary conditions --                  ------------ */

  geometry->x_boundary = PERIODIC;
  geometry->y_boundary = PERIODIC;

  result = xdr_vector(&xdrs, (char *) bc, 2,
		       sizeof(int), (xdrproc_t) xdr_int);
  for (n = 0;  n < 2;  n++)
    geometry->z_boundary_value[n] = (enum boundval) bc[n];

  switch (geometry->z_boundary_value[TOP]) {
  case IRRADIATED:   break;
  case ZERO:         break;
  case THERMALIZED:  break;
  default:
    sprintf(messageStr, "Invalid boundary condition value at TOP: %d",
	    geometry->z_boundary_value[TOP]);
    Error(ERROR_LEVEL_2, routineName, messageStr);
    break;
  }
  switch (geometry->z_boundary_value[BOTTOM]) {
  case IRRADIATED:   break;
  case ZERO:         break;
  case THERMALIZED:  break;
  default:
    sprintf(messageStr, "Invalid boundary condition at BOTTOM: %d",
	    geometry->z_boundary_value[BOTTOM]);
    Error(ERROR_LEVEL_2, routineName, messageStr);
    break;
  }

  /* --- Get horizontal grid spacing and vertical stratification --  */

  result = xdr_double(&xdrs, &geometry->dx);
  geometry->dx *= KM_TO_M;
  result = xdr_double(&xdrs, &geometry->dy);
  geometry->dy *= KM_TO_M;

  geometry->z = (double *) malloc(Nz * sizeof(double));
  result = xdr_vector(&xdrs, (char *) geometry->z, Nz,
		       sizeof(double), (xdrproc_t) xdr_double);
  for (k = 0;  k < Nz;  k++) geometry->z[k] *= KM_TO_M;

  /* --- Allocate space for atmospheric physics --     -------------- */

  atmos->T     = (double *) malloc(Nspace * sizeof(double));
  atmos->ne    = (double *) malloc(Nspace * sizeof(double));
  atmos->vturb = (double *) malloc(Nspace * sizeof(double));

  result = xdr_vector(&xdrs, (char *) atmos->T, Nspace,
		       sizeof(double), (xdrproc_t) xdr_double);
  result = xdr_vector(&xdrs, (char *) atmos->ne, Nspace,
		       sizeof(double), (xdrproc_t) xdr_double);
  result = xdr_vector(&xdrs, (char *) atmos->vturb, Nspace,
		       sizeof(double), (xdrproc_t) xdr_double);

  geometry->vx = (double *) malloc(Nspace * sizeof(double));
  geometry->vy = (double *) malloc(Nspace * sizeof(double));
  geometry->vz = (double *) malloc(Nspace * sizeof(double));

  result = xdr_vector(&xdrs, (char *) geometry->vx, Nspace,
		       sizeof(double), (xdrproc_t) xdr_double);
  result = xdr_vector(&xdrs, (char *) geometry->vy, Nspace,
		       sizeof(double), (xdrproc_t) xdr_double);
  result = xdr_vector(&xdrs, (char *) geometry->vz, Nspace,
		       sizeof(double), (xdrproc_t) xdr_double);

  for (k = 0;  k < Nspace;  k++) {
    atmos->vturb[k] *= KM_TO_M;
    geometry->vx[k] *= KM_TO_M;
    geometry->vy[k] *= KM_TO_M;
    geometry->vz[k] *= KM_TO_M;
  }
  /* --- Determine whether line profiles should account for shift due
         to macroscopic velocities --                    ------------ */

  atmos->moving = FALSE;
  for (k = 0;  k < Nspace;  k++) {
    if ((fabs(geometry->vx[k]) >= atmos->vmacro_tresh) || 
	(fabs(geometry->vy[k]) >= atmos->vmacro_tresh) ||
        (fabs(geometry->vz[k]) >= atmos->vmacro_tresh)) {
      atmos->moving = TRUE;
      break;
    }
  }
  /* --- Finally, get the hydrogen density --          -------------- */

  atmos->nH = matrix_double(atmos->NHydr, Nspace);
  result = xdr_vector(&xdrs, (char *) atmos->nH[0], atmos->NHydr*Nspace,
		      sizeof(double), (xdrproc_t) xdr_double);
  atmos->nHtot = (double *) calloc(Nspace, sizeof(double));
  for (n = 0;  n < atmos->NHydr;  n++) {
    for (k = 0;  k < Nspace;  k++) {
      atmos->nHtot[k] += atmos->nH[n][k];
    }
  }
  /* --- Get angle-quadrature and copy geometry independent quantities
         Nrays and wmu to atmos structure. --          -------------- */

  getAngleQuadr(geometry);
  atmos->Nrays = geometry->Nrays;
  atmos->wmu   = geometry->wmu;

  /* --- Magnetic field is read here --                -------------- */

  atmos->Stokes = readB(atmos);

  getCPU(2, TIME_POLL, "Read Atmosphere");
}
/* ------- end ---------------------------- readAtmos.c ------------- */

/* ------- begin -------------------------- getBoundary.c ----------- */

void getBoundary(Atmosphere *atmos, Geometry *geometry)
{
  const char routineName[] = "getBoundary";

  bool_t result = TRUE;
  int    Nx = geometry->Nx, Ny = geometry->Ny, Nz = geometry->Nz,
         Nplane = geometry->Nplane, Nspect = spectrum.Nspect;

  /* --- Read boundary conditions and irradiation. First VERTICAL
         boundary conditions --                        -------------- */

  switch (geometry->z_boundary_value[TOP]) {
  case IRRADIATED:
    geometry->Itop = matrix_double(Nspect, Nplane);
    result &= xdr_vector(&xdrs, (char *) geometry->Itop[0],
			 Nspect*Nplane, sizeof(double),
			 (xdrproc_t) xdr_double);
    break;
  default:
    geometry->Itop = NULL;
  }

  switch (geometry->z_boundary_value[BOTTOM]) {
  case IRRADIATED:
    geometry->Ibottom = matrix_double(Nspect, Nplane);
    result &= xdr_vector(&xdrs, (char *) geometry->Ibottom[0],
			 Nspect*Nplane, sizeof(double),
			 (xdrproc_t) xdr_double);
    break;
  default:
    geometry->Ibottom = NULL;
  }

  xdr_destroy(&xdrs);
  fclose(atmos->fp_atmos);
}
/* ------- end ---------------------------- getBoundary.c ----------- */
