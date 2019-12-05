/* ------- file: -------------------------- readj.c -----------------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Mon Jan 16 19:56:13 2012 --

       --------------------------                      ----------RH-- */

/* --- Routines to read and write angle-averaged mean intensity and
       background opacities and emissivities. --       -------------- */
 

/* --- Note: pread and pwrite are only able to deal with chunks
             that are less than 2^31 - 1 bytes. The routines
             pread_rh and pwrite_rh read/write chunks of PWRITE_SLICE
             sequentially to overcome this. --         -------------- */

#include <unistd.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "spectrum.h"
#include "background.h"
#include "inputs.h"
#include "error.h"

#define PWRITE_SLICE 1073741824


/* --- Function prototypes --                          -------------- */

size_t pwrite_rh(int fd, const void *buf, size_t count, off_t offset);
size_t pread_rh(int fd, const void *buf, size_t count, off_t offset);


/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern Spectrum spectrum;
extern InputData input;
extern char messageStr[];


/* ------- begin -------------------------- readJlambda.c ----------- */

void readJlambda(int nspect, double *J)
{
  const char routineName[] = "readJlambda";

  bool_t result = TRUE;
  size_t recordsize;
  off_t  offset;

  recordsize = atmos.Nspace * sizeof(double);
  offset     = recordsize * nspect;

  result &= (pread(spectrum.fd_J, J, recordsize, offset) == recordsize);

  if (!result) {
    sprintf(messageStr,
	    "Error reading file: offset = %lld, recordsize = %zu",
	    (long long) offset, recordsize);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }
}
/* ------- end ---------------------------- readJlambda.c ----------- */

/* ------- begin -------------------------- writeJlambda.c ---------- */

void writeJlambda(int nspect, double *J)
{
  const char routineName[] = "writeJlambda";

  bool_t result = TRUE;
  size_t recordsize;
  off_t  offset;

  recordsize = atmos.Nspace * sizeof(double);
  offset     = recordsize * nspect;

  result &= (pwrite(spectrum.fd_J, J, recordsize, offset) == recordsize);

  if (!result) {
    sprintf(messageStr,
	    "Error writing file: offset = %lld, recordsize = %zu",
	    (long long) offset, recordsize);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }
}
/* ------- end ---------------------------- writeJlambda.c ---------- */

/* ------- begin -------------------------- readJ20lambda.c --------- */

void readJ20lambda(int nspect, double *J20)
{
  const char routineName[] = "readJ20lambda";

  bool_t result = TRUE;
  size_t recordsize;
  off_t  offset;

  recordsize = atmos.Nspace * sizeof(double);
  offset     = recordsize * nspect;

  result &= (pread(spectrum.fd_J20, J20,
		   recordsize, offset) == recordsize);

  if (!result) {
    sprintf(messageStr,
	    "Error reading file: offset = %lld, recordsize = %zu",
	    (long long) offset, recordsize);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }
}
/* ------- end ---------------------------- readJ20lambda.c --------- */

/* ------- begin -------------------------- writeJ20lambda.c -------- */

void writeJ20lambda(int nspect, double *J20)
{
  const char routineName[] = "writeJ20lambda";

  bool_t result = TRUE;
  size_t recordsize;
  off_t  offset;

  recordsize = atmos.Nspace * sizeof(double);
  offset     = recordsize * nspect;

  result &= (pwrite(spectrum.fd_J20, J20,
		    recordsize, offset) == recordsize);

  if (!result) {
    sprintf(messageStr,
	    "Error writing file: offset = %lld, recordsize = %zu",
	    (long long) offset, recordsize);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }
}
/* ------- end ---------------------------- writeJ20lambda.c -------- */

/* ------- begin -------------------------- readImu.c --------------- */

void readImu(int nspect, int mu, bool_t to_obs, double *I)
{
  const char routineName[] = "readImu";

  bool_t result = TRUE;
  int    index;
  size_t recordsize;
  off_t  offset;

  recordsize = atmos.Nspace * sizeof(double);

  index  = spectrum.PRDindex[nspect];
  offset = 2*(index*atmos.Nrays + mu) * recordsize;
  if (to_obs)  offset += recordsize;

  result &= (pread(spectrum.fd_Imu, I, 
		    recordsize, offset) == recordsize);

 if (!result) {
    sprintf(messageStr,
	    "Error reading file: offset = %lld, recordsize = %zu",
	    (long long) offset, recordsize);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }
}
/* ------- end ---------------------------- readImu.c --------------- */

/* ------- begin -------------------------- writeImu.c -------------- */

void writeImu(int nspect, int mu, bool_t to_obs, double *I)
{
  const char routineName[] = "writeImu";

  bool_t result = TRUE;
  int    index;
  size_t recordsize;
  off_t  offset;

  recordsize = atmos.Nspace * sizeof(double);

  index  = spectrum.PRDindex[nspect];
  offset = 2*(index*atmos.Nrays + mu) * recordsize;
  if (to_obs)  offset += recordsize;

  result &= (pwrite(spectrum.fd_Imu, I, 
		    recordsize, offset) == recordsize);

 if (!result) {
    sprintf(messageStr,
	    "Error writing file: offset = %lld, recordsize = %zu",
	    (long long) offset, recordsize);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }
}
/* ------- end ---------------------------- writeImu.c -------------- */

/* ------- begin -------------------------- readBackground.c -------- */

void readBackground(int nspect, int mu, bool_t to_obs)
{
  const char routineName[] = "readBackground";

  int    recordno, NrecStokes, NskipStokes;
  size_t recordsize;
  off_t  offset;
  bool_t result = TRUE;
  ActiveSet *as;

  /* --- Read background opacity, emissivity and scattering opacity
         into ActiveSet structure as. This cannot be done sequentially
         with relative offsets because the file may be read by
         different threads at the same time. Therefore, we use absolute
         offsets here and in writeBackground. --       -------------- */

  as = &spectrum.as[nspect];

  recordsize = atmos.Nspace * sizeof(double);

  if (atmos.moving || atmos.Stokes)
    recordno = atmos.backgrrecno[2*(nspect*atmos.Nrays + mu) + to_obs];
  else
    recordno = atmos.backgrrecno[nspect];
  offset = recordno * recordsize;

  /* --- Read emissivity and opacity for Q, U, and V only if we are
         solving explicitly for all four Stokes quantities, but always
         skip the proper amount of records --           ------------- */

  NrecStokes = 1;
  if (atmos.backgrflags[nspect].ispolarized) {
    NskipStokes = 4;
    if (input.StokesMode == FULL_STOKES) NrecStokes = 4;
  } else
    NskipStokes = 1;

  /* --- Read background opacity --                    -------------- */

  result &= (pread_rh(atmos.fd_background, as->chi_c,
		      NrecStokes * recordsize,
		      offset) == NrecStokes * recordsize);
  offset += NskipStokes * recordsize;

  /* --- Read off-diagonal elements propagation matrix K -- --------- */

  if (atmos.backgrflags[nspect].ispolarized && input.magneto_optical) {
    if (input.StokesMode == FULL_STOKES)
      result &= (pread_rh(atmos.fd_background, as->chip_c, 3*recordsize,
			  offset) == 3*recordsize);
    offset += 3 * recordsize;
  }
  /* --- Read background emissivity --                 -------------- */

  result &= (pread_rh(atmos.fd_background, as->eta_c,
                   NrecStokes * recordsize,
                   offset) == NrecStokes * recordsize);
  offset += NskipStokes * recordsize;

  /* --- Read background scattering opacity --         -------------- */

  result &= (pread_rh(atmos.fd_background,
		      as->sca_c, recordsize, offset) == recordsize);

  /* --- Exit if reading is unsuccessful --            -------------- */

  if (!result) Error(ERROR_LEVEL_2, routineName, "Error reading file");
}
/* ------- end ---------------------------- readBackground.c -------- */

/* ------- begin -------------------------- writeBackground.c ------- */

int writeBackground(int nspect, int mu, bool_t to_obs,
                    double *chi_c, double *eta_c, double *sca_c,
                    double *chip_c)
{
  const char routineName[] = "writeBackground";

  int    recordno, Nwrite, NrecStokes;
  bool_t result = TRUE;
  size_t recordsize;
  off_t  offset;

  /* --- Writes background opacity, emissivity, and scattering 
         opacity. Returns Nwrite, the number of records written,
         with a record length of atmos.Nspace * sizeof(double) -- --- */

  recordsize = atmos.Nspace * sizeof(double);

  if (atmos.moving || atmos.Stokes)
    recordno = atmos.backgrrecno[2*(nspect*atmos.Nrays + mu) + to_obs];
  else
    recordno = atmos.backgrrecno[nspect];
  offset = recordno * recordsize;

  if (atmos.backgrflags[nspect].ispolarized)
    NrecStokes = 4;
  else
    NrecStokes = 1;

  Nwrite = 0;

  /* --- Write background opacity --                   -------------- */

  result &= (pwrite_rh(atmos.fd_background, chi_c,
		       NrecStokes * recordsize,
		       offset) == NrecStokes * recordsize);
  Nwrite += NrecStokes;
  offset += NrecStokes * recordsize;

  /* --- Write off-diagonal elements of propagation matrix K -- ----- */

  if (atmos.backgrflags[nspect].ispolarized && input.magneto_optical) {
    result &= (pwrite_rh(atmos.fd_background, chip_c, 3*recordsize,
			 offset) == 3*recordsize);
    Nwrite += 3;
    offset += 3 * recordsize;
  }
  /* --- Write background emissivity --                -------------- */

  result &= (pwrite_rh(atmos.fd_background, eta_c, 
		       NrecStokes * recordsize,
		       offset) == NrecStokes * recordsize);
  Nwrite += NrecStokes;
  offset += NrecStokes * recordsize;

  /* --- Write background scattering opacity --        -------------- */

  result &= (pwrite_rh(atmos.fd_background, sca_c, recordsize,
		       offset) == recordsize);
  Nwrite += 1;

  if (!result)
    Error(ERROR_LEVEL_2, routineName, "Error writing file");

  /* --- Return the number of written records --       -------------- */

  return Nwrite;
}
/* ------- end ---------------------------- writeBackground.c ------- */

/* ------- begin -------------------------- readProfile.c ----------- */

void readProfile(AtomicLine *line, int lamu, double *phi)
{
  const char routineName[] = "readProfile";

  int    Nrecphi, NrecSkip;
  bool_t result = TRUE;
  size_t recordsize;
  off_t  offset;

  if (line->polarizable && (input.StokesMode > FIELD_FREE)) {
    if (input.magneto_optical)
      NrecSkip = 7;
    else
      NrecSkip = 4;
  } else
    NrecSkip = 1;

  if (line->polarizable && input.StokesMode == FULL_STOKES) {
    if (input.magneto_optical)
      Nrecphi = 7;
    else
      Nrecphi = 4;
  } else
    Nrecphi = 1;
  
  recordsize = Nrecphi * atmos.Nspace * sizeof(double);
  offset     = NrecSkip * atmos.Nspace * sizeof(double) * lamu;

  result &= (pread(line->fd_profile, phi, recordsize, offset) ==
	     recordsize);

  if (!result) Error(ERROR_LEVEL_2, routineName, "Error reading file");
}
/* ------- end ---------------------------- readProfile.c ----------- */

/* ------- begin -------------------------- writeProfile.c ---------- */

void writeProfile(AtomicLine *line, int lamu, double *phi)
{
  const char routineName[] = "writeProfile";

  int    Nrecphi;
  bool_t result = TRUE;
  size_t recordsize;
  off_t  offset;

  if (line->polarizable && (input.StokesMode > FIELD_FREE)) {
    if (input.magneto_optical)
      Nrecphi = 7;
    else
      Nrecphi = 4;
  } else
    Nrecphi = 1;

  recordsize = Nrecphi * atmos.Nspace * sizeof(double);
  offset     = recordsize * lamu;

  result &= (pwrite(line->fd_profile, phi, recordsize, offset) ==
	     recordsize);

  if (!result) Error(ERROR_LEVEL_2, routineName, "Error writing file");
}
/* ------- end ---------------------------- writeProfile.c ---------- */

/* ------- begin -------------------------- pread_rh.c -------------- */

size_t pread_rh(int fd, const void *buffer, size_t count, off_t offset)
{
  register int i;

  int    Nslice;
  size_t Nread = 0, remainder;
  off_t  offs = offset;
  char  *buf = (char *) buffer;

  Nslice    = count / PWRITE_SLICE;
  remainder = count % PWRITE_SLICE;

  for (i = 0;  i < Nslice;  i++) {
    Nread += pread(fd, buf, PWRITE_SLICE, offs);
    offs  += PWRITE_SLICE;
    buf   += PWRITE_SLICE;
   }
  if (remainder > 0)
    Nread += pread(fd, buf, remainder, offs);

  return Nread;
}
/* ------- end ---------------------------- pread_rh.c -------------- */

/* ------- begin -------------------------- pwrite_rh.c ------------- */

size_t pwrite_rh(int fd, const void *buffer, size_t count, off_t offset)
{
  register int i;

  int    Nslice;
  size_t Nwrite = 0, remainder;
  off_t  offs = offset;
  char  *buf = (char *) buffer;

  Nslice    = count / PWRITE_SLICE;
  remainder = count % PWRITE_SLICE;

  for (i = 0;  i < Nslice;  i++) {
    Nwrite += pwrite(fd, buf, PWRITE_SLICE, offs);
    offs += PWRITE_SLICE;
    buf  += PWRITE_SLICE;
   }
  if (remainder > 0)
    Nwrite += pwrite(fd, buf, remainder, offs);

  return Nwrite;
}
/* ------- end ---------------------------- pwrite_rh.c ------------- */
