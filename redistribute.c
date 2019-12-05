/* ------- file: -------------------------- redistribute.c ----------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Wed Apr  1 14:01:22 2009 --

       --------------------------                      ----------RH-- */

/* --- Administers iterations to redistribute intensity in PRD line while
       keeping the population number fixed. --         -------------- */
 
#include <stdlib.h>
#include <math.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "spectrum.h"
#include "accelerate.h"
#include "error.h"
#include "inputs.h"


/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern Spectrum spectrum;
extern InputData input;
extern CommandLine commandline;
extern char messageStr[];


/* ------- begin -------------------------- Redistribute.c ---------- */

void Redistribute(int NmaxIter, double iterLimit)
{
  const char routineName[] = "Redistribute";
  register int kr, nact;

  bool_t  quiet, accel, eval_operator, redistribute;
  enum    Interpolation representation;
  int     niter, Nlamu;
  double  drho, drhomax, drhomaxa;
  Atom *atom;
  AtomicLine *line;

  for (nact = 0;  nact < atmos.Nactiveatom;  nact++) {
    atom = atmos.activeatoms[nact];
    
    /* --- Initialize structures for Ng acceleration PRD iteration -- */

    for (kr = 0;  kr < atom->Nline;  kr++) {
      line = &atom->line[kr];
      if (line->PRD && line->Ng_prd == NULL) {
	if (input.PRD_angle_dep)
	  Nlamu = 2*atmos.Nrays * line->Nlambda * atmos.Nspace;
	else
	  Nlamu = line->Nlambda*atmos.Nspace;
	
	line->Ng_prd = NgInit(Nlamu, input.PRD_Ngdelay,
			      input.PRD_Ngorder, input.PRD_Ngperiod,
			      line->rho_prd[0]);
      }
    }
  }
  /* --- Iterate over scattering integral while keeping populations
         fixed --                                      -------------- */

  niter = 1;
  while (niter <= NmaxIter) {

    drhomaxa = 0.0;
    for (nact = 0;  nact < atmos.Nactiveatom;  nact++) {
      atom = atmos.activeatoms[nact];

      drhomax = 0.0;
      for (kr = 0;  kr < atom->Nline;  kr++) {
	line = &atom->line[kr];
	if (line->PRD) {
	  if (input.PRD_angle_dep)
	    PRDAngleScatter(line, representation=LINEAR);
	  else
	    PRDScatter(line, representation=LINEAR);

	  accel = Accelerate(line->Ng_prd, line->rho_prd[0]);
	  sprintf(messageStr, "  PRD: iter #%d, atom %s, line %d,",
		  line->Ng_prd->count-1, atom->ID, kr);
	  drho = MaxChange(line->Ng_prd, messageStr, quiet=FALSE);
	  sprintf(messageStr, (accel) ? " (accelerated)\n" : "\n");
	  Error(MESSAGE, routineName, messageStr);

	  drhomax = MAX(drho, drhomax);
	}
	drhomaxa = MAX(drhomax, drhomaxa);
      }
    }
    /* --- Solve transfer equation with fixed populations -- -------- */

    solveSpectrum(eval_operator=FALSE, redistribute=TRUE);

    if (drhomaxa < iterLimit) break;
    niter++;
  }
}
/* ------- end ---------------------------- Redistribute.c ---------- */
