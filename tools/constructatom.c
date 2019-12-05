/* ------- file: -------------------------- constructatom.c ---------

       Author:        Han Uitenbroek  (huitenbroek@cfa.harvard.edu)
       Last modified: Thu Feb  3 09:40:34 2000 --

       --------------------------                      ----------RH-- */

#include <malloc.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "spectrum.h"
#include "constant.h"
#include "inputs.h"
#include "error.h"
#include "statistics.h"
#include "atomweights.h"

#define COMMENT_CHAR       "#"
#define RLK_RECORD_LENGTH  160
#define RLK_LABEL_SIZE     10
#define Q_WING             20.0
#define MILLI              1.0E-03
#define ANGSTROM_TO_NM     0.1


/* --- The format of the Kurucz linelists are 160 char long:

       FORMAT(F11.4,F7.3,F6.2,F12.3,F5.2,1X,A10,F12.3,F5.2,1X,A10,
              3F6.2,A4,2I2,I3,F6.3,I3,F6.3,2I5,1X,A1,A1,1X,A1,A1,i1,A3.2I5,I6) 
    --                                                 --------------- */

typedef struct {
  char   label[RLK_LABEL_SIZE+1];
  int    lineNo;
  double E, g;
} RLK_level;

int  RLK_ascend(const void *v1, const void *v2);
void writeModelAtom(Atom *atom, FILE *fp_out);


enum Topology topology = ONE_D_PLANE;

Atmosphere atmos;
ProgramStats stats;
InputData input;
CommandLine commandline;
char messageStr[MAX_LINE_SIZE];


/* ------- begin -------------------------- constructatom.c --------- */

void main(int argc, char *argv[])
{
  register int n;

  char   inputLine[RLK_RECORD_LENGTH+1], Gvalues[18+1],
        *commentChar = COMMENT_CHAR, *labeli, *labelj;
  int    Nline, Nread, Nrequired, checkPoint, hfs_i, hfs_j, gL_i, gL_j,
         iso_dl, elemNo, *ll_index, *ul_index, lineNo, NmaxLevel;
  double elem_code, lambda_air, Ji, Jj, Grad, GStark, GvdWaals, pti, Ej,
         lambda, C;
  FILE  *fp_kurucz;
  Atom atom;
  RadTrans *line;
  RLK_Line *rlk_lines, *kptr;
  RLK_level *rlk_levels;

  C = 2 * PI * (Q_ELECTRON/EPSILON_0) * (Q_ELECTRON/M_ELECTRON) / CLIGHT;

  commandline.quiet = FALSE;
  commandline.logfile = stderr;

  if (argc < 2) {
    sprintf(messageStr, "\nUsage: %s filename  NmaxLevel\n", argv[0]);
    message(messageStr);
    exit(0);
  }
  if ((fp_kurucz = fopen(argv[1], "r")) == NULL)
    Error(UNABLE_TO_OPEN_INPUTFILE, argv[0], argv[1]);

  if (argc >= 3)
    NmaxLevel = atof(argv[2]);
  else
    NmaxLevel = 0;

  /* --- Count the number of lines --                  -------------- */

  Nline = 0;
  while (fgets(inputLine, RLK_RECORD_LENGTH+1, fp_kurucz) != NULL)
    if (*inputLine != *commentChar) Nline++;
  rewind(fp_kurucz);

  rlk_lines  = (RLK_Line *)  malloc(Nline*sizeof(RLK_Line));
  rlk_levels = (RLK_level *) malloc(2*Nline * sizeof(RLK_level));

  labeli = (char *) calloc(RLK_LABEL_SIZE+1, sizeof(char));
  labelj = (char *) calloc(RLK_LABEL_SIZE+1, sizeof(char));

  /* --- Read lines from file --                       -------------- */

  kptr   = rlk_lines;
  lineNo = 0;
  while (fgets(inputLine, RLK_RECORD_LENGTH+1, fp_kurucz) != NULL) {
    if (*inputLine != *commentChar) {
      Nread = sscanf(inputLine, "%lf %lf %lf %lf",
		     &lambda_air, &kptr->gf, &elem_code, &kptr->Ei);

      /* --- Verify the element ID --                  -------------- */

      if (kptr == rlk_lines) {
	elemNo = (int) ceil(elem_code);
	sprintf(messageStr, "Found element: %2s\n",
		atomweight[elemNo - 1].ID);
        message(messageStr);
      } else {
	if (elemNo != (int) ceil(elem_code)) {
	  sprintf(messageStr, "Element ID is not the same for all lines.\n"
		  "Found %2s, not %2s\n", atomweight[(int) ceil(elem_code)-1],
		  atomweight[elemNo - 1].ID);
          Error(ERROR_READING_INPUTFILE, argv[0], messageStr);
	}
      }

      air_to_vacuum(1, &lambda_air, &kptr->lambda);
      kptr->gf  = POW10(kptr->gf);
      kptr->stage = (int) (100.0 * modf(elem_code, &pti));
      kptr->pt_index = (int) pti;

      Nread += sscanf(inputLine+36, "%lf", &Ji);
      Nread += sscanf(inputLine+52, "%lf %lf", &Ej, &Jj);
      kptr->gi = 2*Ji + 1;
      kptr->gj = 2*Jj + 1;

      /* --- Convert to Joule --                       -------------- */

      kptr->Ei *= (HPLANCK * CLIGHT) / CM_TO_M;
      Ej *= (HPLANCK * CLIGHT) / CM_TO_M;

      strncpy(labeli, inputLine+42, RLK_LABEL_SIZE);
      strncpy(labelj, inputLine+70, RLK_LABEL_SIZE);

      strncpy(Gvalues, inputLine+80, 18);
      Nread += sscanf(Gvalues, "%lf %lf %lf", &Grad, &GStark, &GvdWaals);
      if (Grad) kptr->Grad = POW10(Grad);
      if (GStark) kptr->GStark = POW10(GStark) * CUBE(CM_TO_M);
      if (GvdWaals) kptr->GvdWaals = POW10(GvdWaals) * CUBE(CM_TO_M);

      Nread += sscanf(inputLine+107, "%d", &kptr->isotope);
      Nread += sscanf(inputLine+109, "%lf", &kptr->hfs_frac);
      kptr->hfs_frac = POW10(kptr->hfs_frac);
      Nread += sscanf(inputLine+118, "%lf", &kptr->iso_frac);
      kptr->iso_frac = POW10(kptr->iso_frac);
      Nread += sscanf(inputLine+124, "%5d%5d", &hfs_i, &hfs_j);
      kptr->hfs_i = ((double) hfs_i) * MILLI * KBOLTZMANN;
      kptr->hfs_j = ((double) hfs_j) * MILLI * KBOLTZMANN;

      Nread += sscanf(inputLine+143, "%5d%5d", &gL_i, &gL_j);
      kptr->gL_i = gL_i * MILLI;
      kptr->gL_j = gL_j * MILLI;

      /* Nread += sscanf(inputLine+155, "%d", &iso_dl); */
      iso_dl = 0;
      kptr->iso_dl = iso_dl * MILLI * ANGSTROM_TO_NM;

      rlk_levels[2*lineNo].lineNo = lineNo;
      rlk_levels[2*lineNo].E = fabs(kptr->Ei);
      rlk_levels[2*lineNo].g = kptr->gi;
      strcpy(rlk_levels[2*lineNo].label, labeli); 

      rlk_levels[2*lineNo + 1].lineNo = (lineNo) ? -lineNo : -Nline;
      rlk_levels[2*lineNo + 1].E = fabs(Ej);
      rlk_levels[2*lineNo + 1].g = kptr->gj;
      strcpy(rlk_levels[2*lineNo + 1].label, labelj); 

      checkNread(Nread, Nrequired=17, argv[0], checkPoint=1);
      kptr++;
      lineNo++;
    }
  }
  fclose(fp_kurucz);

  strcpy(atom.ID, atomweight[elemNo-1].ID);
  atom.weight = atomweight[elemNo-1].weight;
  sprintf(messageStr, "Read %d Kurucz lines for element %2s from file %s\n",
	  Nline, atom.ID, argv[1]);
  message(messageStr);

  message("Sorting energy levels...\n");
  qsort((void *) rlk_levels, 2*Nline, sizeof(struct RLK_level),
	RLK_ascend);

  ll_index = (int *) calloc(Nline, sizeof(int));
  ul_index = (int *) calloc(Nline, sizeof(int));

  /* --- Collect unique levels as distinguished by energy and
         statistical weight (in case of degeneracy) -- -------------- */

  message("Collecting unique levels... \n");
  atom.Nlevel = 1;

  for (n = 1;  n < 2*Nline;  n++) {
    if (rlk_levels[n].g != rlk_levels[atom.Nlevel-1].g || 
	strcmp(rlk_levels[n].label, rlk_levels[atom.Nlevel-1].label)) {

      rlk_levels[atom.Nlevel].E = rlk_levels[n].E;
      rlk_levels[atom.Nlevel].g = rlk_levels[n].g;
      rlk_levels[atom.Nlevel].lineNo = rlk_levels[n].lineNo;
      strcpy(rlk_levels[atom.Nlevel].label, rlk_levels[n].label);

      atom.Nlevel++;
    }
    if (rlk_levels[n].lineNo < 0)
      ul_index[(-rlk_levels[n].lineNo) % Nline] = atom.Nlevel-1;
    else
      ll_index[rlk_levels[n].lineNo] = atom.Nlevel-1;
  }
  sprintf(messageStr, "Found %d unique levels\n", atom.Nlevel);
  message(messageStr);

  if (NmaxLevel > 0) atom.Nlevel = MIN(atom.Nlevel, NmaxLevel);
  atom.E = (double *) malloc(atom.Nlevel * sizeof(double));
  atom.g = (double *) malloc(atom.Nlevel * sizeof(double));
  atom.stage = (int *)   malloc(atom.Nlevel * sizeof(int));
  atom.label = (char **) malloc(atom.Nlevel * sizeof(char *));

  for (n = 0;  n < atom.Nlevel;  n++) {
    atom.E[n] = rlk_levels[n].E;
    atom.g[n] = rlk_levels[n].g;

    if (rlk_levels[n].lineNo < 0)
      atom.stage[n] = rlk_lines[(-rlk_levels[n].lineNo) % Nline].stage;
    else
      atom.stage[n] = rlk_lines[rlk_levels[n].lineNo].stage;

    atom.label[n] = (char *) calloc(ATOM_LABEL_WIDTH+1, sizeof(char));
    switch (atom.stage[n]) {
    case 0: sprintf(atom.label[n], "SI I ");    break;
    case 1: sprintf(atom.label[n], "SI II ");   break;
    case 2: sprintf(atom.label[n], "SI III ");  break;
    case 3: sprintf(atom.label[n], "SI IV ");   break;
    case 4: sprintf(atom.label[n], "SI V ");    break;
    case 5: sprintf(atom.label[n], "SI VI ");   break;
    }
    strcat(atom.label[n], rlk_levels[n].label);
  }
  atom.Nfixed = atom.Ncont = atom.Nprd = 0;

  atom.radtrans = (struct RadTrans *) malloc(Nline * sizeof(struct RadTrans));
  line = atom.radtrans;
  atom.Nline = 0;
  kptr = rlk_lines;
  for (n = 0;  n < Nline;  n++) {
    initRadTrans(line);
    line->type = BOUNDBOUND;
    line->parent.atom = &atom;
    line->lambda0 = kptr->lambda;

    line->i = ll_index[n];
    line->j = ul_index[n];

    if (line->i < atom.Nlevel  &&  line->j < atom.Nlevel) {
      lambda = kptr->lambda * NM_TO_M;
      line->Aji = C / SQ(lambda) * kptr->gf / kptr->gj;
      line->Bji = CUBE(lambda) / (2.0 * HPLANCK * CLIGHT) * line->Aji;
      line->Bij = (kptr->gj / kptr->gi) * line->Bji;

      line->shape = VOIGT;
      line->Nlambda = 15;
      line->qcore = 1.0;
      line->qwing = 20.0;
      line->lvdW = UNSOLD;
      line->cvdWaals[2] = line->cvdWaals[0] = 1.0;
      line->cvdWaals[3] = line->cvdWaals[1] = 0.0;
      line->cStark = 1.0;
      line->Grad = kptr->Grad;

      line->symmetric = TRUE;
      atom.Nline++;
      line++;
    }
    kptr++;
  }
  atom.radtrans = (struct RadTrans *)
    realloc(atom.radtrans, atom.Nline * sizeof(struct RadTrans));

  message("Writing model atom....\n");
  writeModelAtom(&atom, stdout);
}
/* ------- end ---------------------------- constructatom.c --------- */

/* ------- begin -------------------------- RLK_ascend.c ------------ */

int RLK_ascend(const void *v1, const void *v2)
{
  struct RLK_level *s1 = (struct RLK_level *) v1,
                   *s2 = (struct RLK_level *) v2;

  if (s1->E < s2->E)
    return -1;
  else if (s1->E > s2->E)
    return 1;
  else {

    /* --- In case of energy degeneracy order in ascending g -- ----- */

    if (s1->g < s2->g)
      return -1;
    else if (s1->g > s2->g)
      return 1;
    else
      return 0;
  }
}
/* ------- end ---------------------------- RLK_ascend.c ------------ */
