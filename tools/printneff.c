/* ------- file: -------------------------- printneff.c ------------- */

/* Reads atomic model and prints labels and corresponding effective
 * quantum numbers.
 *
 * Han Uitenbroek
 * Last modified: Mar 6, 1996
 */
 
#include <math.h>
#include <stdlib.h>

#include "rh.h"
#include "atom.h"
#include "constant.h"
#include "error.h"
#include "inputs.h"

/* --- Function prototypes --                          -------------- */

void rawAtom(Atom *atom, char *atomFileName);


/* --- Global variables --                             -------------- */

CommandLine commandline;
char   messageStr[MAX_LINE_SIZE];


/* ------- begin -------------------------- printneff.c ------------- */

int main( int argc, char *argv[] )
{
  register int i;

  char   errorMsg[MAX_LINE_SIZE];
  int    ic;
  float  n_eff, Z, E_Rydberg;
  Atom atom;
  FILE  *fpOut;

  commandline.quiet = FALSE;
  commandline.logfile = stderr;

  if (argc >= 3)
    fpOut = fopen( argv[2], "w" );
  else if (argc == 2)
    fpOut = stdout;
  else {
    fprintf( stderr, "Usage: %s inFile [outFile]\n", argv[0] );
    exit( 0 );
  }
  rawAtom(&atom, argv[1]);

  /* --- Rydberg constant for finite atomic mass --    -------------- */

  E_Rydberg = E_RYDBERG / (1.0 + M_ELECTRON / (atom.weight * AMU));

  fprintf(fpOut, "\n       label              n_eff\n");
  fprintf(fpOut, "--------------------------------\n");
  for (i = 0;  i < atom.Nlevel-1;  i++) {
    for (ic = i + 1;
	 ((atom.stage[ic] < atom.stage[i]+1) && (ic < atom.Nlevel));  ic++);
    if (atom.stage[ic] == atom.stage[i]) {
      sprintf(errorMsg, "Found no overlying continuum for level %d", i);
      Error(ERROR_LEVEL_2, argv[0], errorMsg);
    } else {
      Z = (float) (atom.stage[i] + 1);
      n_eff = Z * sqrt(E_Rydberg / (atom.E[ic] - atom.E[i]));
      fprintf(fpOut, "'%20s'  %f\n", atom.label[i], n_eff);
    }
  }
}
/* ------- end ---------------------------- printneff.c ------------- */
