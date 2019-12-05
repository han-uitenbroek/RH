/* ------- file: -------------------------- duplicate.c -------------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Tue Jan 25 16:17:47 2000 --

       --------------------------                      ----------RH-- */

/* --- Routine duplicateLevel checks whether the multiplet level
       designated by label is present in the atomic model atom.

       Routine duplicateLine checks whether line transition between
       levels designated by labeli and labelj is present in the
       atomic model atom.

 Note: Ambiguity due to fine-splitting is prevented by cutting the label
       off at the term designation (ie. by comparing up till the odd or
       even designation). This requires consistent naming of the levels
       in active and background models of the same species of course!
       --                                              -------------- */
 
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "background.h"
#include "error.h"

/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

extern char messageStr[];


/* ------- begin -------------------------- duplicateLevel.c -------- */

bool_t duplicateLevel(Atom *atom, char *label)
{
  const char routineName[] = "duplicateLevel";
  register int i;

  char   multiplet[ATOM_LABEL_WIDTH+1], *ptr;
  bool_t duplicate;
  int    length;

  if (atom == NULL) return FALSE;

  strcpy(multiplet, label);
  ptr = multiplet + (strlen(multiplet) - 1);
  while ((*ptr != 'E')  &&  (*ptr != 'O')  &&  (ptr > multiplet))  ptr--;
  if (ptr > multiplet)
    *(ptr + 1) = '\0';
  else {
    sprintf(messageStr, "Cannot determine parity of atomic level %s", label);
    Error(ERROR_LEVEL_1, routineName, messageStr);
    return FALSE;
  }

  length = strlen(multiplet);
  duplicate = FALSE;
  for (i = 0;  i < atom->Nlevel;  i++) {

    /* --- Compare up to the length of the multiplet -- ------------- */

    if (strncmp(atom->label[i], multiplet, length) == 0) {
      duplicate = TRUE;
      break;
    }
  }
  return duplicate;
}
/* ------- end ---------------------------- duplicateLevel.c -------- */

/* ------- begin -------------------------- duplicateLine.c --------- */

bool_t duplicateLine(struct Atom *atom, char *labeli, char *labelj)
{
  const char routineName[] = "duplicateLine";
  register int kr;

  char   multipleti[ATOM_LABEL_WIDTH+1], multipletj[ATOM_LABEL_WIDTH+1],
        *ptr;
  bool_t duplicate;
  int    i, j, lengthi, lengthj;

  if (atom == NULL) return FALSE;

  strcpy(multipleti, labeli);
  ptr = multipleti + (strlen(multipleti) - 1);
  while ((*ptr != 'E')  &&  (*ptr != 'O')  &&  (ptr > multipleti))  ptr--;
  if (ptr > multipleti)
    *(ptr + 1) = '\0';
  else {
    sprintf(messageStr,
	    "Cannot determine parity of atomic level %s", labeli);
    Error(ERROR_LEVEL_1, routineName, messageStr);
  }
  lengthi = strlen(multipleti);

  strcpy(multipletj, labelj);
  ptr = multipletj + (strlen(multipletj) - 1);
  while ((*ptr != 'E')  &&  (*ptr != 'O')  &&  (ptr > multipletj))  ptr--;
  if (ptr > multipletj)
    *(ptr + 1) = '\0';
  else {
    sprintf(messageStr,
	    "Cannot determine parity of atomic level %s", labelj);
    Error(ERROR_LEVEL_1, routineName, messageStr);
  }
  lengthj = strlen(multipletj);

  duplicate = FALSE;
  for (kr = 0;  kr < atom->Nline;  kr++) {
    i = atom->line[kr].i;
    j = atom->line[kr].j;
    if (strncmp(atom->label[i], multipleti, lengthi) == 0  &&
	strncmp(atom->label[j], multipletj, lengthj) == 0) {
      duplicate = TRUE;
      break;
    }
  }
  return duplicate;
}
/* ------- end ---------------------------- duplicateLine.c --------- */
