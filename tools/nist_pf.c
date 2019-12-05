#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "constant.h"
#include "inputs.h"
#include "error.h"

#define MAX_REC_LENGTH  142
#define COMMENT_CHAR    "#"

struct NIST_Level {
  char  label[ATOM_LABEL_WIDTH+1];
  int   parity, stage;
  float E, g, g_Landee;
};

CommandLine commandline;
char   messageStr[MAX_LINE_SIZE];

int main(int argc, char *argv[])
{
  register int i, k;

  char   inputLine[MAX_REC_LENGTH+1], *commentChar = COMMENT_CHAR,
         *atomID, *spectrum, *label, *ptr, *Jstr, *term1, *term2,
        *Estr;
  int    Nlevel, Nread, Ntemp, n, J;
  float  E, g, *kT, *pf;
  FILE  *fp_NIST;
  struct NIST_Level *levels, *nl;

  commandline.quiet = FALSE;
  commandline.logfile = stderr;

  if (argc < 2) {
    fprintf(stderr, "Usage: %s inFile E0\n", argv[0]);
    exit(0);
  }
  /* --- Open the data file for model atom --          -------------- */

  if ((fp_NIST = fopen(argv[1], "r")) == NULL) {
    sprintf(messageStr, "Unable to open input file %s", argv[1]);
    Error(ERROR_LEVEL_2, argv[0], messageStr);
  }

  Nlevel = 0;
  while (fgets(inputLine, MAX_REC_LENGTH, fp_NIST) != NULL)
    if (*inputLine != *commentChar) Nlevel++;
  printf("Found %d levels in NIST input file %s\n\n", Nlevel, argv[1]);
  rewind(fp_NIST);

  levels = (struct NIST_Level *) malloc(Nlevel * sizeof(struct NIST_Level));

  nl = levels;
  while (fgets(inputLine, MAX_REC_LENGTH+1, fp_NIST) != NULL &&
	  nl < levels + Nlevel-1) {
    if (*inputLine != *commentChar) {
      atomID   = strtok(inputLine, " ");
      spectrum = strtok(NULL, " ");
      label    = strtok(NULL, " ");
      UpperCase(atomID);
      for (i = 0;  i < strlen(label);  i++)
	if (*(label+i) == '.') *(label+i) = ' ';
      UpperCase(label);

      if (strcmp(spectrum, "I") == 0)
	nl->stage = 0;
      else if (strcmp(spectrum, "II") == 0)
        nl->stage = 1;
      else if (strcmp(spectrum, "III") == 0)
        nl->stage = 2;

      term1 = strtok(NULL, " ");
      term2 = strtok(NULL, " ");

      Jstr  = strtok(NULL, " ");
      if (strstr(Jstr, "/")) {
	sscanf(Jstr, "%d/", &J);
        nl->g = J + 1.0;
      } else {
	sscanf(Jstr, "%d", &J);
	nl->g = 2.0*J + 1.0;
      }
      /* --- Parity, @ if odd, blank if even --        -------------- */

      Estr = strtok(NULL, " ");
      if (strstr(Estr, "@")) {
	nl->parity = 1;
        Estr = strtok(NULL, " ");
      } else
	nl->parity = 0;
      nl->E = atof(Estr) * (HPLANCK * CLIGHT) / CM_TO_M;

      sprintf(nl->label, "%2s %s %-8.8s %2.2s%1s",
	      atomID, spectrum, label, term2, (nl->parity) ? "O" : "E");
      nl++;
    }
  }
  Ntemp = argc - 2;

  kT = (float *) malloc(Ntemp * sizeof(float));
  for (k = 0;  k < Ntemp;  k++)
    kT[k] = 1.0 / (KBOLTZMANN * atof(argv[k+2]));
  pf = (float *) calloc(Ntemp, sizeof(float));

  for (i = 0;  i < Nlevel-1;  i++) {
    nl = levels + i;
    if (nl->stage == levels[0].stage) {
      for (k = 0;  k < Ntemp;  k++) {
	/* printf("%s -> %e\n", nl->label, nl->g * exp(-nl->E * kT[k])); */
	pf[k] += nl->g * exp(-nl->E * kT[k]);
      }
    }
  }

  for (k = 0;  k < Ntemp;  k++)
    printf(" T = %f     pf = %e\n", 1.0/(KBOLTZMANN*kT[k]), pf[k]);
}
