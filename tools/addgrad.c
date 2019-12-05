/* ------- file: -------------------------- addgrad.c --------------- */

/* adds Aji for appropriate transitions to get natural broadening.
 *
 * Han Uitenbroek
 * Last modified: Thu Feb  3 09:41:49 2000 --
 */
 
#include <stdlib.h>

#include "rh.h"
#include "atom.h"
#include "error.h"
#include "inputs.h"

/* --- Function prototypes --                          -------------- */

void rawAtom(Atom *atom, char *atomFileName);
void writeModelAtom(Atom *atom, FILE *fp_out);


/* --- Global variables --                             -------------- */

CommandLine commandline;
char messageStr[MAX_LINE_SIZE];


/* ------- begin -------------------------- addgrad.c --------------- */

int main( int argc, char *argv[] )
{
  FILE  *fp_out;
  Atom atom;

  commandline.quiet = FALSE;
  commandline.logfile = stderr;

  if (argc >= 3)
    fp_out = fopen( argv[2], "w" );
  else if (argc == 2)
    fp_out = stdout;
  else {
    fprintf( stderr, "Usage: %s inFile [outFile]\n", argv[0] );
    exit( 0 );
  }
  rawAtom(&atom, argv[1]);
  writeModelAtom(&atom, fp_out);
}
/* ------- end ---------------------------- addgrad.c --------------- */
