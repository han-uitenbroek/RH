/* ------- file: -------------------------- waveinfo.c --------------

       Version:       rh1.0
       Author:        Han Uitenbroek (huitenbroek@sunspot.noao.edu)
       Last modified: Mon Jul  8 14:28:37 2002 --

       --------------------------                      ----------RH-- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "spectrum.h"
#include "xdr.h"


#define CHECK_OPTION(string, option, lmin) \
        (strncmp((string), (option), MAX((lmin),strlen(string))) == 0)

struct Option {
  char *name;
  int   deflt;
  int   minlen;
  char *cmnt;
};

/* ------- begin -------------------------- ShowUsage.c ------------- */

void ShowUsage(char *programName, struct Option *theOpt)
{
  fprintf(stderr, " Usage:  %s [Options] Filename\n", programName);
  fprintf(stderr, "  where Options can be any of the following:\n\n");

  while (theOpt->name != NULL) {
    fprintf(stderr, "   %s \t<%s> (default: %d)\n", 
	    theOpt->name, theOpt->cmnt, theOpt->deflt);
    theOpt++;
  }
  fprintf(stderr, 
  "\n  Each option can be abreviated to its unambiguous length.\n\n");

  exit(-1);
}
/* ------- end ---------------------------- ShowUsage.c ------------- */

/* ------- begin -------------------------- waveinfo.c -------------- */

int main(int argc, char *argv[])
{
  register int i;
  const char routineName[] = "waveinfo";

  char    format[] = "%9.3f   %9.3f   %9.3f   %9.3f   %9.3f\n",
           messageStr[132];
  bool_t  showfull = FALSE, to_air = FALSE, result = TRUE;
  int     Nlambda;
  double *lambda;
  FILE   *fp_wave;
  XDR     xdrs;

  /* --- Reads wavelength input file and diplays contents -- -------- */

  struct Option theOpt[] = {
    {"-help",  FALSE, 2, "Print this message          "},
    {"-full",  FALSE, 2, "Print all wavelengths       "},
    {"-air",   FALSE, 2, "Convert to wavelength in air"},
    NULL
  };
  if (argc <= 1) ShowUsage(argv[0], theOpt);

  for (i = 1;  i < argc  &&  argv[i][0] == '-';  i++) {

    if(CHECK_OPTION(argv[i], theOpt[0].name, theOpt[0].minlen)) {
      ShowUsage(argv[0], theOpt);
    }
    else if(CHECK_OPTION(argv[i], theOpt[1].name, theOpt[1].minlen)) {
      showfull = TRUE;
    }
    else if(CHECK_OPTION(argv[i], theOpt[2].name, theOpt[2].minlen)) {
      to_air = TRUE;
    }
    else {
      fprintf(stderr, " Invalid or ambiguous option: %s\n", argv[i]); 
      ShowUsage(argv[0], theOpt);
    }
  }

  if ((fp_wave = fopen(argv[argc-1], "r")) == NULL) {
    fprintf(stderr, " Error opening file: %s\n", argv[argc-1]);
    exit(0);
  }
  xdrstdio_create(&xdrs, fp_wave, XDR_DECODE);

  result &= xdr_int(&xdrs, &Nlambda);

  lambda = (double *) malloc(Nlambda * sizeof(double));
  result &= xdr_vector(&xdrs, (char *) lambda, Nlambda,
		       sizeof(double), (xdrproc_t) xdr_double);

  if (!result) {
    sprintf(messageStr, "Unable to read from input file %s",
	    argv[argc-1]);
  }

  xdr_destroy(&xdrs);
  fclose(fp_wave);

  if (to_air) vacuum_to_air(Nlambda, lambda, lambda);

  fprintf(stderr,
	  "\n Found %d wavelengths, ranging from %9.3f to %9.3f [nm]%s\n\n",
	  Nlambda, lambda[0], lambda[Nlambda-1],
	  (to_air) ? " (in air)" : " (in vaccuum)");
  if (showfull) {
    for (i = 0;  i < Nlambda-5;  i += 5) {
      printf(format, lambda[i], lambda[i+1], lambda[i+2],
	     lambda[i+3], lambda[i+4]);
    }
    for (i = MAX(Nlambda-5, 0);  i < Nlambda;  i++)
      printf("%9.3f%s", lambda[i], (i == Nlambda-1) ? "\n\n" : "   ");
  }
  free(lambda);
  exit(1);
}
/* ------- end ---------------------------- waveinfo.c -------------- */
