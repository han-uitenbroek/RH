/* ------- file: -------------------------- ludcmp.c ----------------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Wed Apr 19 15:23:29 2000 --

       --------------------------                      ----------RH-- */

/* --- Matrix inversion based on LU decomposition, from Press, Flannery,
       Teukolsky and Vetterling, in Numerical Recipes.

 Note: The routines here are modified so that they adhere to the more
       natural mapping of A_ij = A[i][j], which is the transposed
       of the convention used by Press et al. --       -------------- */


#include <stdlib.h>
#include <math.h>

#include "rh.h"
#include "error.h"

/* --- Function prototypes --                          -------------- */

void LUdecomp(int N, double **A, int *index, double *d);
void LUbacksubst(int N, double **A, int *index, double *b);


/* --- Global variables --                             -------------- */

extern char messageStr[];


/* ------- begin -------------------------- SolveLinearEq.c --------- */

void SolveLinearEq(int N, double **A, double *b, bool_t improve)
{
  register int i, j;

  int   *index;
  double d, **A_copy, *b_copy, *residual;

  /* --- If improve == TRUE improve the solution of the set of
         linear equations by evaluating the residual and correcting
	 the initial solution.

    See: Press, Flannery, Teukolsky and Vetterling, Numerical Recipes,
         The art of scientific computing 1986, p. 41
          --                                           -------------- */

  index = (int *) malloc(N * sizeof(int));

  /* --- Copy matrix and source vector --              -------------- */

  if (improve) {
    residual = (double *) malloc(N * sizeof(double));
    b_copy = (double *) malloc(N * sizeof(double));
    A_copy = matrix_double(N, N);

    for (i = 0;  i < N;  i++) {
      b_copy[i] = b[i];
      for (j = 0;  j < N;  j++) A_copy[i][j] = A[i][j];
    }
  }
  /* --- Initial solution --                           ------------- */

  LUdecomp(N, A, index, &d);
  LUbacksubst(N, A, index, b);

  if (improve) {
    for (i = 0;  i < N;  i++) {
      residual[i] = b_copy[i];
      for (j = 0;  j < N;  j++) residual[i] -= A_copy[i][j] * b[j];
    }
    LUbacksubst(N, A, index, residual);

    /* --- Correct the initial solution --             ------------- */

    for (i = 0;  i < N;  i++) b[i] += residual[i];

    free(residual);
    free(b_copy);
    freeMatrix((void **) A_copy);
  }
  free(index);
}

/* ------- end ---------------------------- SolveLinearEq.c --------- */

/* ------- begin -------------------------- LUdecomp.c -------------- */

#define TINY 1.0e-20;

void LUdecomp(int N, double **A, int *index, double *d)
{
  register int i, j, k;

  int imax = 0;
  double big, dum, sum, temp, *vv;

  vv = (double *) malloc(N * sizeof(double));
  *d = 1.0;

  for (i = 0;  i < N;  i++) {
    big = 0.0;
    for (j = 0;  j < N;  j++)
      if ((temp = fabs(A[i][j])) > big)  big = temp;
    if (big == 0.0) {
      sprintf(messageStr, "Singular matrix");
      Error(ERROR_LEVEL_2, "LUdecomp", messageStr);
    }
    vv[i] = 1.0 / big;
  }

  for (j = 0;  j < N;  j++) {
    for (i = 0;  i < j;  i++) {
      sum = A[i][j];
      for (k = 0;  k < i;  k++) sum -= A[i][k] * A[k][j];
      A[i][j] = sum;
    }
    big = 0.0;
    for (i = j;  i < N;  i++) {
      sum = A[i][j];
      for (k = 0;  k < j;  k++)
	sum -= A[i][k] * A[k][j];
      A[i][j] = sum;
      if ((dum = vv[i]*fabs(sum)) >= big) {
	big  = dum;
	imax = i;
      }
    }
    if (j != imax) {
      for (k = 0;  k < N;  k++) {
	dum = A[imax][k];
	A[imax][k] = A[j][k];
	A[j][k] = dum;
      }
      *d = -(*d);
      vv[imax] = vv[j];
    }
    index[j] = imax;
    if (A[j][j] == 0.0) A[j][j] = TINY;
    if (j != N) {
      dum = 1.0 / A[j][j];
      for (i = j+1;  i < N;  i++) A[i][j] *= dum;
    }
  }
  free(vv);
}

/* ------- end ---------------------------- LUdecomp.c -------------- */

/* ------- begin -------------------------- LUbacksubst.c ----------- */

void LUbacksubst(int N, double **A, int *index, double *b)
{
  register int i, j;

  int    ii = -1, ip;
  double sum;
 
  for (i = 0;  i < N;  i++) {
    ip = index[i];
    sum = b[ip];
    b[ip] = b[i];
    if (ii >= 0) {
      for (j = ii;  j < i;  j++) sum -= A[i][j] * b[j];
    } else if (sum)
      ii = i;
    b[i] = sum;
  }
  for (i = N-1;  i >= 0;  i--) {
    sum = b[i];
    for (j = i+1;  j < N;  j++) sum -= A[i][j]*b[j];
    b[i] = sum / A[i][i];
  }
}
/* ------- end ---------------------------- LUbacksubst.c ----------- */
