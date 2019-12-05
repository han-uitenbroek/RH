/* ------- file: -------------------------- accelerate.c ------------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Mon Apr 17 12:00:40 2000 --

       --------------------------                      ----------RH-- */

/* --- Implements Ng convergence acceleration with general order.

  See: K. C. Ng 1974, J. Chem. Phys. 61, 2680
 Also: L. Auer 1987, in "Numerical Radiative Transfer",
         ed. W. Kalkofen, pp. 101-109
       --                                              -------------- */

#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#include "rh.h"
#include "accelerate.h"
#include "error.h"
#include "statistics.h"


/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */


/* ------- begin -------------------------- NgInit.c ---------------- */

struct Ng* NgInit(int N, int Ndelay, int Norder, int Nperiod,
		  double *solution)
{
  /* --- Initialize data structure and allocate space for previous
         solutions, coefficient matrix A and correction vector b,
         copy initial solution into previous solution matrix -- ----- */

  register int k;

  struct Ng *Ngs;

  Ngs = (struct Ng*) malloc(sizeof(struct Ng));

  Ngs->N       = N;
  Ngs->Norder  = Norder;
  Ngs->Nperiod = Nperiod;
  Ngs->Ndelay  = MAX(Ndelay, Norder + 2);

  if (Norder > 0) {
    Ngs->A = matrix_double(Norder, Norder);
    Ngs->b = (double *) malloc(Norder * sizeof(double));
  }
  Ngs->previous   = matrix_double(Norder + 2, N);
  Ngs->theStorage = Ngs->previous[0];
  for (k = 0;  k < N;  k++)
    Ngs->previous[0][k] = solution[k];
  Ngs->count = 1;

  return Ngs;
}
/* ------- end ---------------------------- NgInit.c ---------------- */

/* ------- begin -------------------------- Accelerate.c ------------ */

bool_t Accelerate(struct Ng *Ngs, double *solution)
{
  register int i, j, k;

  int      Norder = Ngs->Norder, ip, ipp, i0;
  double **Delta, *weight;

  /* --- Store the current solution --                --------------- */

  i = Ngs->count % (Norder + 2);
  for (k = 0;  k < Ngs->N;  k++) {
    Ngs->previous[i][k] = solution[k];
  }
  (Ngs->count)++;

  /* --- Accelerate only after we have accumulated enough iterations,
         and then only every Ngs->period th iteration after iteration
         Ngs->Ndelay --                               --------------- */

  if ((Norder > 0) && (Ngs->count >= Ngs->Ndelay)  &&
      !((Ngs->count - Ngs->Ndelay) % Ngs->Nperiod)) {
    getCPU(4, TIME_START, NULL);

    Delta  = matrix_double(Norder+1, Ngs->N);
    weight = (double *)  malloc(Ngs->N * sizeof(double));

    for (i = 0;  i <= Norder;  i++) {
      ip  = (Ngs->count - 1 - i) % (Norder + 2);
      ipp = (Ngs->count - 2 - i) % (Norder + 2);
      for (k = 0;  k < Ngs->N;  k++) {
	Delta[i][k] = Ngs->previous[ip][k] - Ngs->previous[ipp][k];
      }
    }
    /* --- Use weighted acceleration --                -------------- */

    for (k = 0;  k < Ngs->N;  k++)  weight[k] = 1.0 / fabs(solution[k]);
    for (i = 0;  i < Norder;  i++) {
      Ngs->b[i] = 0.0;
      for (j = 0;  j < Norder;  j++)  Ngs->A[i][j] = 0.0;
    }
    
    /* --- Fill the coefficient matrix and invert linear set --  ---- */

    for (j = 0;  j < Norder;  j++) {
      for (k = 0;  k < Ngs->N;  k++) {
	Ngs->b[j] += weight[k] * Delta[0][k]*(Delta[0][k] - Delta[j+1][k]);
      }
      for (i = 0;  i < Norder;  i++) {
	for (k = 0;  k < Ngs->N;  k++) {
	  Ngs->A[i][j] += weight[k] *
	    (Delta[j+1][k] - Delta[0][k]) * (Delta[i+1][k] - Delta[0][k]);
	}
      }
    }
    SolveLinearEq(Norder, Ngs->A, Ngs->b, TRUE);

      /* --- Construct the linear combination for the accelerated
	     solution, and restore the accelerated solution --  ----- */

    i0 = (Ngs->count - 1) % (Norder + 2);
    for (i = 0;  i < Norder;  i++) {
      ip = (Ngs->count - 2 - i) % (Norder + 2);
      for (k = 0;  k < Ngs->N;  k++) {
	solution[k] += Ngs->b[i] *
	  (Ngs->previous[ip][k] - Ngs->previous[i0][k]);
      }
    }
    for (k = 0;  k < Ngs->N;  k++) Ngs->previous[i0][k] = solution[k];

    /* --- Clean up the temporary variables --        --------------- */

    free(weight);
    freeMatrix((void **) Delta);

    getCPU(4, TIME_POLL, "Accelerate");
    return TRUE;
  } else
    return FALSE;
}
/* ------- end ---------------------------- Accelerate.c ------------ */

/* ------- begin -------------------------- NgFree.c ---------------- */

void NgFree(struct Ng *Ngs)
{
  freeMatrix((void **) Ngs->previous);

  if (Ngs->Norder > 0) {
    free(Ngs->b);
    freeMatrix((void **) Ngs->A);
  }

  free(Ngs);
}
/* ------- end ---------------------------- NgFree.c ---------------- */
