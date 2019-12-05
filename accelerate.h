/* ------- begin -------------------------- Accelerate.h ------------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Tue Dec 28 14:02:56 1999 --

       --------------------------                      ----------RH-- */

#ifndef __ACCELERATE_H__
#define __ACCELERATE_H__


/* --- Defines structure and prototypes for Ng acceleration -- ------ */

struct Ng {
  int      N, Ndelay, Norder, Nperiod, count;
  double **previous, **A, *b, *theStorage;
};


/* --- Associated function prototypes --               -------------- */

bool_t Accelerate(struct Ng *Ngs, double *solution);
void   NgFree(struct Ng *Ngs);
struct Ng *NgInit(int N, int Ndelay, int Norder, int Nperiod,
		  double *solution);
double MaxChange(struct Ng *Ngs, char *text, bool_t quiet);


#endif /* !__ACCELERATE_H__ */

/* ------- end ---------------------------- Accelerate.h ------------ */
