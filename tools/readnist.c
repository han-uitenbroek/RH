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
char messageStr[MAX_LINE_SIZE];

int main(int argc, char *argv[])
{
  register int i;

  char   inputLine[MAX_REC_LENGTH+1], *commentChar = COMMENT_CHAR,
         atomID[ATOM_ID_WIDTH+1], spectrum[6], label[25], *ptr, term[3],
         Jstr[5];
  int    Nlevel, Nread, n, J;
  float  E0, E, g;
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
  E0 = atof(argv[2]);

  Nlevel = 0;
  while (fgets(inputLine, MAX_REC_LENGTH, fp_NIST) != NULL)
    if (*inputLine != *commentChar) Nlevel++;
  printf("Found %d levels in NIST input file %s\n\n", Nlevel, argv[1]);
  rewind(fp_NIST);

  levels = (struct NIST_Level *) malloc(Nlevel * sizeof(struct NIST_Level));

  nl = levels;
  while (fgets(inputLine, MAX_REC_LENGTH+1, fp_NIST) != NULL) {
    if (*inputLine != *commentChar) {
      Nread = sscanf(inputLine, "%2s %s %s",
		     atomID, spectrum, label);
      UpperCase(atomID);
      for (i = 0;  i < strlen(label);  i++)
	if (*(label+i) == '.') *(label+i) = ' ';
      UpperCase(label);
      sscanf(inputLine+37, "%2s", term);

      if (strcmp(spectrum, "I") == 0)
	nl->stage = 0;
      else if (strcmp(spectrum, "II") == 0)
        nl->stage = 1;
      else if (strcmp(spectrum, "III") == 0)
        nl->stage = 2;

      /* --- Parity, @ if odd, blank if even --        -------------- */

      nl->parity = (*(inputLine+57) == '@') ? 1 : 0;
      sprintf(nl->label, "%2s %s %-9.9s %2s%1s", atomID, spectrum, label,
	      term, (nl->parity) ? "O" : "E");

      Nread += sscanf(inputLine+41, "%s", Jstr);
      if (strstr(Jstr, "/")) {
	sscanf(Jstr, "%d/", &J);
        nl->g = J + 1.0;
      } else {
	sscanf(Jstr, "%d", &J);
	nl->g = 2.0*J + 1.0;
      }
      Nread += sscanf(inputLine+60, "%f %f", &nl->E, &nl->g_Landee);
      nl++;
    }
  }
  fprintf(stdout,
	  "%1s   E[cm^-1]    g           label[20]         stage   levelNo\n",
          COMMENT_CHAR);
  fprintf(stdout, "%1s                     '|----|----|----|----'\n",
          COMMENT_CHAR);

  strncpy(label, levels[0].label, ATOM_LABEL_WIDTH);
  g = E = 0.0;
  n = 0;
  for (i = 0;  i < Nlevel-1;  i++) {
    g += levels[i].g;
    E += levels[i].E * levels[i].g;
    if (strncmp(levels[i+1].label, label, ATOM_LABEL_WIDTH) != 0) {
      fprintf(stdout, " %10.3f  %6.2f   '%-20s'    %1d     %3d\n",
	    E0 + E/g, g, levels[i].label, levels[i].stage, n);
      g = E = 0.0;
      n++;
      strncpy(label, levels[i+1].label, ATOM_LABEL_WIDTH);
    }
  }
  fprintf(stdout, " %10.3f  %6.2f   '%-20s'    %1d     %3d\n",
	  E0 + levels[Nlevel-1].E, levels[Nlevel-1].g,
	  levels[Nlevel-1].label, levels[Nlevel-1].stage, n);
}
