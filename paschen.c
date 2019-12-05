/* ------- file: -------------------------- paschen.c ---------------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Wed Aug  5 09:57:21 2009 --

       --------------------------                      ----------RH-- */

  /* --- Routines for implementing the incomplete Paschen-Back effect.

    See: E. Landi Degl'Innocenti & M. Landolfi, 2004
              "Polarization in Spectral Lines", Kluwer

         Section 3.4
         --                                            -------------- */

#include <stdlib.h>
#include <math.h>

#include "rh.h"
#include "atom.h"
#include "error.h"
#include "constant.h"


#define N_MAX_ITER_TQLI 30

#define SIGN(a, b) ((b) >= 0.0 ? fabs(a) : -fabs(a))


/* --- Function prototypes --                          -------------- */

double Pythag(double a, double b);
void   tqli(double *d, double *e, int N, double **z);
void   Paschen(double L, double S, double *E0, double B,
	       Paschenstruct *ps);


/* --- Global variables --                             -------------- */

extern char messageStr[];

/* ------- begin -------------------------- Paschen_Back.c ---------- */

void Paschen_Back()
{

}
/* ------- end ---------------------------- Paschen_Back.c ---------- */

/* ------- begin -------------------------- Paschen.c --------------- */

void Paschen(double L, double S, double *E0, double B, Paschenstruct *ps)
{
  register int M, j, jp;

  int Mindex;
  double Larmor, *ud, J, Jp, Jmin, Jmax, Ja, Jb;

  /* --- See: E. Landi Degl'Innocenti & M. Landolfi, 2004
              "Polarization in Spectral Lines", Kluwer

         Eqs. 3.57 and 3.58
         --                                            -------------- */

  /* --- Larmor frequency in energy units (J) --       -------------- */

  Larmor = (HPLANCK * Q_ELECTRON / (4.0 * PI *M_ELECTRON)) * B;

  Jmin = fabs(L - S);
  Jmax = L + S;

  Mindex = 0;
  for (M = -Jmax;  M <=Jmax;  M++) {
    Ja = MAX(fabs(M), Jmin);
    Jb = Jmax;
    ps[Mindex].Nj = (int)(Jb - Ja) + 1;

    ps[Mindex].eigenval = (double *) malloc(ps[Mindex].Nj * sizeof(double));
    ps[Mindex].C = matrix_double(ps[Mindex].Nj, ps[Mindex].Nj);

    ud = (double *) malloc(ps[Mindex].Nj * sizeof(double));

    for (j = 0;  j < ps[Mindex].Nj;  j++) {
      J = Ja + j;
 
      ps[Mindex].eigenval[j] += E0[j] +
	Larmor * (M + (((int) (2*J + L + S + M)) % 2 ? -1.0 : 1.0) *
		  (2.0*J + 1.0) * sqrt(S * (S + 1.0) * (2.0*S + 1.0)) *
		  w3js(J, J, 1.0, -M, M, 0.0) * w6js(J, J, 1.0, S, S, L));
    }

    if (ps[Mindex].Nj > 1) {
      for (jp = 1;  jp < ps[Mindex].Nj;  jp++) {
	Jp = Ja + jp;
	J  = Jp - 1.0;

	ud[jp] = Larmor * (((int) (J + Jp + L + S + M)) % 2 ? -1.0 : 1.0) *
	  sqrt((2.0*J +1.0) * (2.0*Jp + 1.0) * S * (S + 1.0) * (2.0*S + 1.0)) *
	  w3js(J, Jp, 2.0, -M, M, 0.0) * w6js(J, Jp, 1.0, S, S, L);
      }
    }

    for (j = 0;  j < ps[Mindex].Nj;  j++) ps[Mindex].C[j][j] = 1.0;
    tqli(ps[Mindex].eigenval, ud, ps[Mindex].Nj, ps[Mindex].C);

    free(ud);
    Mindex++;
  }
}
/* ------- end ---------------------------- Paschen.c --------------- */

/* ------- begin -------------------------- Pythag.c ---------------- */

double Pythag(double a, double b)
{
  double absa, absb;

  /* --- Returns sqrt(a^2 + b^2) without destructive
         under- or overflow.

    See: Press, Flannery, Teukolsky and Vetterling, 1992
         Numerical Recipes in C,
         The art of scientific computing, p. 70
          --                                           -------------- */

  absa = fabs(a);
  absb = fabs(b);

  if (absa > absb)
    return absa * sqrt(1.0 + SQ(absb / absa));
  else
    return ((absb == 0.0) ? 0.0 : absb * sqrt(1.0 + SQ(absa / absb)));
}
/* ------- end ---------------------------- Pythag.c ---------------- */

/* ------- begin -------------------------- tqli.c ------------------ */

void tqli(double *d, double *e, int N, double **z)
{
  const char routineName[] = "tqli";
  register int m, l, iter, i, k;

  double s, r, p, g, f, dd, c, b;

  /* --- Tridiagonal QL implicit eigenvector decomposition (tqli.c).

    See: Press, Flannery, Teukolsky and Vetterling, 1992
         Numerical Recipes in C,
         The art of scientific computing, p. 480

    Note: Modified for proper C indexing starting at 0.
          --                                           -------------- */

  if (N == 1) return;

  for (i = 1;  i < N;  i++) e[i - 1] = e[i];
  e[N-1] = 0.0;

  for (l = 0;  l < N;  l++) {
    iter = 0;
    do {
      for (m = l;  m < N - 1;  m++) {
	dd = fabs(d[m]) + fabs(d[m + 1]);
	if ((double)(fabs(e[m]) + dd) == dd) break;
      }

      if (m != l) {
	if (iter++ == N_MAX_ITER_TQLI) {
	  sprintf(messageStr,
		  "  Too many iterations in tqli: %d\n", iter-1);
	  Error(ERROR_LEVEL_2, routineName, messageStr);
	}

	g = (d[l + 1] - d[l]) / (2.0 * e[l]);
	r = Pythag(g, 1.0);
	g = d[m] - d[l] + e[l] / (g + SIGN(r, g));
	c = 1.0;
	s = 1.0; 
	p = 0.0;
	for (i = m-1;  i >= l;  i--) {
	  f = s * e[i];
	  b = c * e[i];
	  r = Pythag(f, g);
	  e[i + 1] = r;

	  if (r == 0.0) {
	    d[i + 1] -= p;
	    e[m] = 0.0;
	    break;
	  }
	  s = f / r;
	  c = g / r;
	  g = d[i + 1] - p;
	  r = (d[i] - g) * s + 2.0 * c * b;
          p = s * r;
	  d[i + 1] = g + p;
	  g = c * r - b;
	  for (k = 0;  k < N;  k++) {
	    f = z[k][i + 1];
	    z[k][i + 1] = s * z[k][i] + c * f;
	    z[k][i]     = c * z[k][i] - s * f;
	  }
	}
	if (r == 0.0  &&  i >= l) continue;
	d[l] -= p;
	e[l]  = g;
	e[m]  = 0.0;
      }
    } while (m != l);
  }
}
/* ------- end ---------------------------- tqli.c ------------------ */
