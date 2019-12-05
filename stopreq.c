/* ------- file: -------------------------- Stoprequested.c ---------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Mon Oct  9 11:27:42 2000 --

       --------------------------                      ----------RH-- */

/* --- Checks whether iterations should be stopped prematurely while
       still allowing intermediate results to be written to output
       in the proper way.

       Current version checks for presence of file STOP_RH in running
       directory as suggested by Jo Bruls (KIS) --     -------------- */


#include <stdio.h>

#include "rh.h"
#include "error.h"

#define STOP_REQUEST_FILE "STOP_RH"


/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

extern char messageStr[];


/* ------- begin -------------------------- StopRequested.c --------- */

bool_t StopRequested()
{
  const char routineName[] = "StopRequested";

  if (fopen(STOP_REQUEST_FILE, "r")) {
    sprintf(messageStr, "Stopping iterations because file %s is present",
	  STOP_REQUEST_FILE);
    Error(WARNING, routineName, messageStr);
    return TRUE;
  } else {
    return FALSE;
  }
}
/* ------- end ---------------------------- StopRequested.c --------- */
