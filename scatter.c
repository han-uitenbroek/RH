/* ------- file: -------------------------- scatter.c ---------------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Thu Feb 22 17:02:07 2024 --

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
#include <unistd.h>
#include <errno.h>
#include <sys/types.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "spectrum.h"
#include "inputs.h"
#include "constant.h"
#include "error.h"
#include "statistics.h"
#include "parallel.h"

#define TENSION  8.0
#define PRD_FILE_TEMPLATE "PRD_%s_%d-%d.dat"


/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern Spectrum spectrum;
extern InputData input;
extern char messageStr[];


/* --- The loop over spatial locations in PRDScatter is parallelized:
       the scattering integral at one location uses only that location's
       data and writes only rho_prd[.][k] for that k. The redistribution
       weights live in a scratch file with one block per location (all
       emission wavelengths of the line, in order), for each subordinate
       line in turn. Each location's block has a known offset, so
       workers read and write their own blocks with pread/pwrite instead
       of sharing a file position, a window of consecutive locations at
       a time: the blocks follow each other in the file, and large
       requests are what the disk and the OS cache handle well. --- */

#define PRD_IO_WINDOW  (2*1024*1024)

typedef struct {
  Atom       *atom;
  AtomicLine *PRDline, **XRD;
  int         Nsubordinate, fd;
  enum Interpolation representation;
  bool_t      initialize;
  double     *Pj, *adamp;
  long       *blocksize;      /* doubles in each block, [kxrd*Nspace + k] */
  off_t      *kstart;         /* byte offset of each block                */
} PRDjob;

static void prdIntegrate(int kbeg, int kend, int nt, void *arg);

/* ------- begin -------------------------- prdGrid.c --------------- */

static int prdGrid(double q_emit, double waveratio, double *q0, double *qN)
{
  /* --- Integration limits over absorption wavelength, using only
	 regions where the redistribution function is non-zero, and the
	 number of points of the fine grid. (See also function GII.) - */

  if (fabs(q_emit) < PRD_QCORE) {
    *q0 = -PRD_QWING;
    *qN =  PRD_QWING;
  } else {
    if (fabs(q_emit) < PRD_QWING) {
      if (q_emit > 0.0) {
	*q0 = -PRD_QWING;
	*qN = waveratio * (q_emit + PRD_QSPREAD);
      } else {
	*q0 = waveratio * (q_emit - PRD_QSPREAD);
	*qN = PRD_QWING;
      }
    } else {
      *q0 = waveratio * (q_emit - PRD_QSPREAD);
      *qN = waveratio * (q_emit + PRD_QSPREAD);
    }
  }
  return (int) ((*qN - *q0) / PRD_DQ) + 1;
}
/* ------- end ---------------------------- prdGrid.c --------------- */

/* ------- begin -------------------------- prdBlockSizes.c --------- */

static void prdBlockSizes(int kbeg, int kend, int nt, void *arg)
{
  PRDjob     *job = (PRDjob *) arg;
  AtomicLine *PRDline = job->PRDline;
  int         kxrd, la, k;
  long        n;
  double      q_emit, q0, qN, waveratio;

  for (kxrd = 0;  kxrd < job->Nsubordinate;  kxrd++) {
    waveratio = job->XRD[kxrd]->lambda0 / PRDline->lambda0;
    for (k = kbeg;  k < kend;  k++) {
      n = 0;
      for (la = 0;  la < PRDline->Nlambda;  la++) {
	q_emit = (PRDline->lambda[la] - PRDline->lambda0) * CLIGHT /
	  (PRDline->lambda0 * job->atom->vbroad[k]);
	n += prdGrid(q_emit, waveratio, &q0, &qN);
      }
      job->blocksize[kxrd*atmos.Nspace + k] = n;
    }
  }
}
/* ------- end ---------------------------- prdBlockSizes.c --------- */

/* ------- begin -------------------------- prdWindow.c ------------- */

/* --- Window of locations k to ke-1, to be read in one go (or filled and
       then written in one go): as many locations as fit in
       PRD_IO_WINDOW bytes, but at least one. --       -------------- */

static int prdWindowEnd(long *blocksize, int k, int kend)
{
  size_t bytes = 0;
  int    ke = k;

  while (ke < kend  &&
	 (ke == k  ||  bytes + blocksize[ke]*sizeof(double) <= PRD_IO_WINDOW)) {
    bytes += blocksize[ke] * sizeof(double);
    ke++;
  }
  return ke;
}

static size_t prdWindowSize(long *blocksize, off_t *kstart, int k, int ke)
{
  return (kstart[ke-1] - kstart[k]) / sizeof(double) + blocksize[ke-1];
}

static void prdWindowIO(PRDjob *job, AtomicLine *XRDline, double *win,
			size_t Nwin, off_t offset, int k0, int k1,
			bool_t write_it)
{
  const char routineName[] = "PRDScatter";

  char   *p = (char *) win;
  size_t  left = Nwin * sizeof(double);
  ssize_t n;

  while (left > 0) {

    /* --- A single call may not take more than INT_MAX bytes on some
	   systems, and may be interrupted --          -------------- */

    size_t chunk = (left < (size_t) 1 << 30) ? left : (size_t) 1 << 30;

    n = (write_it) ? pwrite(job->fd, p, chunk, offset) :
		     pread(job->fd, p, chunk, offset);
    if (n < 0  &&  errno == EINTR) continue;
    if (n <= 0) {
      sprintf(messageStr,
	      "Unable to %s the redistribution weights\n"
	      " Line %d -> %d, locations %d to %d",
	      (write_it) ? "write" : "read", XRDline->j, XRDline->i,
	      k0, k1 - 1);
      Error(ERROR_LEVEL_2, routineName, messageStr);
    }
    p += n;  offset += n;  left -= n;
  }
}
/* ------- end ---------------------------- prdWindow.c ------------- */


/* ------- begin -------------------------- PRDScatter.c ------------ */

void PRDScatter(AtomicLine *PRDline, enum Interpolation representation)
{
  const char routineName[] = "PRDScatter";
  register int  la, k, kr, ip, kxrd;

  char    filename[MAX_LINE_SIZE];
  bool_t  initialize;
  int     ij, Nsubordinate, nthreads;
  off_t   offset;
  double *adamp, cDop, *Pj;
  Atom *atom;
  AtomicLine *line, **XRD;
  AtomicContinuum *continuum;
  PRDjob  job;

  /* --- Routine for angle-averaged and hybrid angle-dependent
         approximation. In the latter Jgas is used instead of J -- -- */ 

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
  /* --- The spline interpolators keep their coefficients in static
	 storage, so only the linear representation runs in parallel - */

  nthreads = (representation == LINEAR) ? input.Nthreads : 1;

  job.atom           = atom;
  job.PRDline        = PRDline;
  job.XRD            = XRD;
  job.Nsubordinate   = Nsubordinate;
  job.representation = representation;
  job.initialize     = initialize;
  job.fd             = fileno(PRDline->fp_GII);
  job.Pj             = Pj;
  job.adamp          = adamp;
  job.blocksize = (long *)  malloc(Nsubordinate*atmos.Nspace * sizeof(long));
  job.kstart    = (off_t *) malloc(Nsubordinate*atmos.Nspace * sizeof(off_t));

  /* --- Where each location's block of weights sits in the file: the
	 size of every block first, then their running sum --  -------- */

  parallelFor(atmos.Nspace, nthreads, prdBlockSizes, &job);
  offset = 0;
  for (kxrd = 0;  kxrd < Nsubordinate;  kxrd++) {
    for (k = 0;  k < atmos.Nspace;  k++) {
      job.kstart[kxrd*atmos.Nspace + k] = offset;
      offset += job.blocksize[kxrd*atmos.Nspace + k] * sizeof(double);
    }
  }

  /* --- Loop over subordinate lines and all spatial locations, with
	 the locations split between the threads --    -------------- */

  parallelFor(atmos.Nspace, nthreads, prdIntegrate, &job);

  /* --- Clean temporary variable space --             -------------- */

  free(job.blocksize);  free(job.kstart);
  free(XRD);    free(Pj);     free(adamp);

  sprintf(messageStr, "Scatter Int %5.1f", PRDline->lambda0);
  getCPU(3, TIME_POLL, messageStr);
}
/* ------- end ---------------------------- PRDScatter.c ------------ */

/* ------- begin -------------------------- prdIntegrate.c ---------- */

static void prdIntegrate(int kbeg, int kend, int nt, void *arg)
{
  register int  la, k, lap, kxrd;

  PRDjob *job = (PRDjob *) arg;
  bool_t  hunt, initialize = job->initialize;
  int     Np, Nsubordinate = job->Nsubordinate, kwin_start, kwin_end;
  long    pos, *blocksize;
  off_t  *kstart;
  size_t  Nwin = 0, wincap = 0;
  enum Interpolation representation = job->representation;
  double  q_emit, q0, qN, *q_abs = NULL, *qp = NULL, *wq = NULL,
	 *qpp = NULL, *gii, *adamp = job->adamp, gnorm, *J = NULL,
	  Jbar, scatInt, *J_k = NULL, *Pj = job->Pj, gamma, waveratio,
	 *win = NULL, *wbuf;
  Atom *atom = job->atom;
  AtomicLine *PRDline = job->PRDline, **XRD = job->XRD, *XRDline;

  /* --- Loop over subordinate lines --                -------------- */

  for (kxrd = 0;  kxrd < Nsubordinate;  kxrd++) {
    XRDline = XRD[kxrd];
    waveratio = XRDline->lambda0 / PRDline->lambda0;

    J_k   = realloc(J_k, XRDline->Nlambda * sizeof(double));
    q_abs = realloc(q_abs, XRDline->Nlambda * sizeof(double));

    /* --- This line's blocks of weights in the file --  -------------- */

    blocksize  = job->blocksize + kxrd*atmos.Nspace;
    kstart     = job->kstart + kxrd*atmos.Nspace;
    kwin_start = kwin_end = kbeg;

    /* --- Loop over the spatial locations of this worker -- --------- */

    for (k = kbeg;  k < kend;  k++) {

      /* --- Start a new window of weights when the previous one is
	     used up (first writing out the finished one when the weights
	     are being computed for the first time) --   -------------- */

      if (k == kwin_end) {
	if (initialize  &&  k > kbeg)
	  prdWindowIO(job, XRDline, win, Nwin, kstart[kwin_start],
		      kwin_start, kwin_end, TRUE);
	kwin_start = k;
	kwin_end   = prdWindowEnd(blocksize, k, kend);
	Nwin       = prdWindowSize(blocksize, kstart, kwin_start, kwin_end);
	if (Nwin > wincap) {
	  wincap = (Nwin > PRD_IO_WINDOW/sizeof(double)) ?
	    Nwin : PRD_IO_WINDOW/sizeof(double);
	  win = (double *) realloc(win, wincap * sizeof(double));
	}
	if (!initialize)
	  prdWindowIO(job, XRDline, win, Nwin, kstart[kwin_start],
		      kwin_start, kwin_end, FALSE);
      }
      wbuf = win + (kstart[k] - kstart[kwin_start]) / sizeof(double);
      pos  = 0;

      gamma = atom->n[XRDline->i][k] / atom->n[PRDline->j][k] *
	XRDline->Bij / Pj[k];
      Jbar = XRDline->Rij[k] / XRDline->Bij;

      /* --- Get local mean intensity and wavelength in Doppler units */

      if (input.PRD_angle_dep == PRD_ANGLE_APPROX) {
	for (la = 0;  la < XRDline->Nlambda;  la++) {
	  J_k[la] = spectrum.Jgas[XRDline->Nblue + la][k];
	  q_abs[la] = (XRDline->lambda[la] - XRDline->lambda0) * CLIGHT /
	    (XRDline->lambda0 * atom->vbroad[k]);
	}
      } else if (input.PRD_angle_dep == PRD_ANGLE_AVER) {
	for (la = 0;  la < XRDline->Nlambda;  la++) {
	  J_k[la] = spectrum.J[XRDline->Nblue + la][k];
	  q_abs[la] = (XRDline->lambda[la] - XRDline->lambda0) * CLIGHT /
	    (XRDline->lambda0 * atom->vbroad[k]);
	}
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
 
	Np = prdGrid(q_emit, waveratio, &q0, &qN);
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

	gii = wbuf + pos;
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
	pos += Np;
      }
    }

    /* --- The last window of weights --               -------------- */

    if (initialize  &&  kend > kbeg)
      prdWindowIO(job, XRDline, win, Nwin, kstart[kwin_start],
		  kwin_start, kend, TRUE);
  }

  free(J_k);    free(q_abs);  free(qp);   free(qpp);
  free(wq);     free(J);      free(win);
}
/* ------- end ---------------------------- prdIntegrate.c ---------- */

/* ------- begin -------------------------- PRDAngleScatter.c ------- */

void PRDAngleScatter(AtomicLine *PRDline,
		     enum Interpolation representation)
{
  const char routineName[] = "PRDAngleScatter";
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
	    (PRDline->polarizable && input.StokesMode == FULL_STOKES))
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
