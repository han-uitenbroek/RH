/* ------- file: -------------------------- hunt.c ------------------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Tue Feb 16 14:40:59 1999 --

       --------------------------                      ----------RH-- */

/* --- Find index of value in array (cf., Num. Recipes, p 91).
       Use previous value of ilow to hunt up or down list and bracket value.
       --                                              -------------- */
 
#include "rh.h"

/* ------- begin -------------------------- Hunt.c ------------------ */

void Hunt(int n, double *array, double value, int *ilow)
{
  bool_t ascend;
  int    ihigh, index, increment;

  ascend = (array[n-1] > array[0]) ? TRUE : FALSE;
  if ((*ilow <= 0)  ||  (*ilow > n-1)) {

    /* --- Input guess not useful here, go to bisection --  --------- */

    *ilow = 0;
    ihigh = n;
  } else {

    /* --- Else hunt up or down to bracket value --    -------------- */ 

    increment = 1;
    if (((value >= array[*ilow]) ? TRUE : FALSE) == ascend) {
      ihigh = *ilow + increment;
      if (*ilow == n-1) return;

      /* --- Hunt up --                                -------------- */

      while (((value >= array[ihigh]) ? TRUE : FALSE) == ascend) {
	*ilow = ihigh;
	increment += increment;
	ihigh = *ilow + increment;
        if (ihigh >= n) { ihigh = n;  break; }
      }
    } else {
      ihigh = *ilow;
      if (*ilow == 0) return;

      /* --- Hunt down --                              -------------- */

      while (((value <= array[*ilow]) ? TRUE : FALSE) == ascend) {
	ihigh = *ilow;
	increment += increment;
	*ilow = ihigh - increment;
        if (*ilow <= 0) { *ilow = 0;  break; }
      }
    }
  }
  /* --- Bisection algorithm --                        -------------- */

  if (ascend) {
    while (ihigh - *ilow > 1) {
      index = (ihigh + *ilow) >> 1;
      if (value >= array[index])
	*ilow = index;
      else
	ihigh = index;
    }
  } else {
    while (ihigh - *ilow > 1) {
      index = (ihigh + *ilow) >> 1;
      if (value <= array[index])
	*ilow = index;
      else
	ihigh = index;
    }
  }
}
/* ------- end ---------------------------- Hunt.c ------------------ */

/* ---------------------------------------- Locate.c ---------------- */

/* --- Find index of value in array (cf., Num. Recipes, p 90).

 Note: The Num. Recipes routine does not give the correct index
       for values that are exactly equal to an array value!
       --                                              -------------- */
 
/* ------- begin -------------------------- Locate.c ---------------- */

void Locate(int n, double *array, double value, int *ilow)
{
  bool_t ascend;
  int    ihigh, index;

  ascend = (array[n-1] > array[0]) ? TRUE : FALSE;
  *ilow = 0;  ihigh = n;

  if (ascend) {
    while (ihigh - *ilow > 1) {
      index = (ihigh + *ilow) >> 1;
      if (value >= array[index])
	*ilow = index;
      else
	ihigh = index;
    }
  } else {
    while (ihigh - *ilow > 1) {
      index = (ihigh + *ilow) >> 1;
      if (value <= array[index])
	*ilow = index;
      else
	ihigh = index;
    }
  }
}
/* ------- end ---------------------------- Locate.c ---------------- */
