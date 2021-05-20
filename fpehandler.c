/* ------- file: -------------------------- fpehandler.c ------------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Thu May 20 11:29:23 2021 --

       --------------------------                      ----------RH-- */

/* --- Trap floating point exceptions on various machines:

       --                                              -------------- */

#include <math.h>
#include <stdio.h>

#include "rh.h"
#include "error.h"

extern char messageStr[];


#if defined(Darwin)

void SetFPEtraps(void)
{
  /* --- Explicitly do not set traps --                -------------- */

  Error(MESSAGE, "SetFPEtraps", 
	"\n Darwin: FPE traps have not been set explicitly\n");
}
/* --- end SETNOTRAPS --                               -------------- */

#elif defined(Linux)


#define _GNU_SOURCE 1
#include <fenv.h>

int feenableexcept(int excepts);

void SetFPEtraps(void)
{
  int result;
  
  /* --- Enable some exceptions.
         At startup all exceptions are masked. --      -------------- */

  result = feenableexcept(FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
}

#else

/* --- unknown --                                      -------------- */

void SetFPEtraps(void)
{
  Error(MESSAGE, "SetFPEtraps",
	"\nUnsupported CPU and/or OS: cannot set FPE traps explicitly\n");
}

#endif

/* ------- end ------------------------------------------------------ */
