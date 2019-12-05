/* ------- file: -------------------------- lande.c ----------------- */

/* Print effective Lande factors for lines in atom
 *
 * Han Uitenbroek
 * Last modified: Thu Feb  3 11:51:26 2000 --
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


/* ------- begin -------------------------- lande.c ----------------- */

void main(int argc, char *argv[])
{
  register int kr;
  Atom atom;
  AtomicLine *line;

  commandline.quiet = FALSE;
  commandline.logfile = stderr;

  if (argc < 2) {
    fprintf(stderr, "Usage: %s atomfile\n", argv[0]);
    exit(0);
  }
  rawAtom(&atom, argv[1]);

  for (kr = 0;  kr < atom.Nline;  kr++) {
    line = atom.line + kr;
    if (line->g_Lande_eff == 0.0) {
      printf("Line %d -> %d:  g_eff = %f\n", line->j, line->i,
	     effectiveLande(line));
    } else {
      printf("Line %d -> %d:  g_eff = %f (set in atom file)\n",
	     line->j, line->i, line->g_Lande_eff);
    }
  }
}
/* ------- end ---------------------------- lande.c ----------------- */
