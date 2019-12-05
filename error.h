/* ---------------------------------------- error.h -----------------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Wed Dec  8 14:57:23 1999 --

       --------------------------                      ----------RH-- */

#ifndef __ERROR_H__
#define __ERROR_H__

/* --- Defines error numbers, messages and associated actions -- ---- */


/* ------- begin -------------------------- error.h ----------------- */

enum errorlevel {MESSAGE, WARNING, ERROR_LEVEL_1, ERROR_LEVEL_2};


/* --- Associated function prototypes --               -------------- */

void Error(enum errorlevel level, const char *routineName,
	   const char *messageStr);

#endif /* !__ERROR_H__ */

/* ------- end ---------------------------- error.h ----------------- */
