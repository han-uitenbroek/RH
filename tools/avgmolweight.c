/* ------- file: -------------------------- avgmolweight.c ----------

       Version:       rh1.0, tools
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Tue Apr 24 17:27:26 2007 --

       --------------------------                      ----------RH-- */

/* --- Reads standard format abundance file and prints average
       molecular weight --                             -------------- */

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "inputs.h"

/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

CommandLine commandline;
InputData input;
char messageStr[MAX_LINE_SIZE];


/* ------- begin -------------------------- avgmolweight.c ---------- */

int main(int argc, char *argv[])
{
  Atmosphere atmos;

  commandline.quiet = FALSE;
  commandline.logfile = stderr;

  if (argc < 2) {
    fprintf(stderr, "Usage: %s abund_file Kurucz_pf_file\n", argv[0]);
    exit(0);
  }
  strcpy(input.abund_input, argv[1]);
  strcpy(input.KuruczData, "none");
  strcpy(input.pfData, argv[2]);
  readAbundance(&atmos);
  fprintf(stderr, "\n  The average weight per Hydrogen atom for \n"
	  "  abundance file %s is %6.4f", argv[1], atmos.wght_per_H);
  fprintf(stderr, "\n  The average molecular weight for \n"
	  "  abundance file %s is %6.4f\n\n", argv[1], atmos.avgMolWght);
}
/* ------- end ---------------------------- avgmolweight.c ---------- */
