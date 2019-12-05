/* ------- file: -------------------------- fpehandler.c ------------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Wed Jan  4 17:47:35 2017 --

       --------------------------                      ----------RH-- */

/* --- Trap floating point exceptions on various machines:

       --                                              -------------- */

#include <math.h>
#include <stdio.h>

#include "rh.h"
#include "error.h"

extern char messageStr[];

#if defined(SETNOTRAPS)

/* --- If SETNOTRAPS has been defined (see Makefile) -- ------------- */

void SetFPEtraps(void)
{
  /* --- Explicitly do not set traps --                -------------- */

  Error(MESSAGE, "SetFPEtraps", 
	"\nFPE traps have not been set explicitly\n");
}
/* --- end SETNOTRAPS --                               -------------- */

#else

/* --- SunOS --                                        -------------- */

#if defined(SunOS)

/* --- Traps floating point exceptions for SunOS5.

   OS: SunOS5

  see: ``man ieee_handler''

       Requires linking with -lsunmath
       --                                              -------------- */

#include <stdlib.h>
#include <sunmath.h>
#include <siginfo.h>
#include <ucontext.h>

void Trapped_FPE_Exception(int sig, siginfo_t *sip, ucontext_t *uap);

void SetFPEtraps(void)
{
  const char routineName[] = "SetFPEtraps";

  if (ieee_handler("set", "common", 
                   (sigfpe_handler_type) Trapped_FPE_Exception) != 0)
    Error(MESSAGE, routineName, "IEEE trapping not supported here\n");
  else
    Error(MESSAGE, routineName,
	  "\n-Setting FPE traps for sparc (SunOS 5.x)\n");
}

void Trapped_FPE_Exception(int sig, siginfo_t *sip, ucontext_t *uap)
{
  const char routineName[] = "Trapped_FPE_Exception";

  char *type = "unknown";
 
  switch(sip->si_code) {
 
  case FPE_INTDIV:   type = "integer divide by zero          "; break;
  case FPE_INTOVF:   type = "integer overflow                "; break;
  case FPE_FLTDIV:   type = "floating point division by zero "; break;
  case FPE_FLTUND:   type = "floating point underflow        "; break;
  case FPE_FLTOVF:   type = "floating point overflow         "; break;
  case FPE_FLTRES:   type = "floating point inexact          "; break;
  case FPE_FLTINV:   type = "invalid floating point operation"; break;
  case FPE_FLTSUB:   type = "subscript out of range          "; break;
  }
 
  sprintf(messageStr, "  ---- trapped IEEE FPE: %s ----\n"
	  "  ---- signal: %d, code: 0x%x, aborting ----\n",
	  type, sig, sip->si_code);
  Error(MESSAGE, routineName, messageStr);
  abort();
}

/* --- end SunOS --                                    -------------- */

/* --- Linux --                                        -------------- */

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

/* --- end Linux --                                    -------------- */

#else

/* --- unknown --                                      -------------- */

void SetFPEtraps(void)
{
  Error(MESSAGE, "SetFPEtraps",
	"\nUnsupported CPU and/or OS: cannot set FPE traps explicitly\n");
}

#endif

/* --- end else --                                     -------------- */

#endif
/* ------- end ------------------------------------------------------ */
