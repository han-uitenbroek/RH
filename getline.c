/* ------- file: -------------------------- getline.c ---------------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Thu Dec  9 11:43:18 1999 --

       --------------------------                      ----------RH-- */

/* --- Routines for reading formatted input files. --  -------------- */

 
#include <ctype.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>

#include "rh.h"
#include "error.h"

/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

extern char messageStr[];


/* ------- begin -------------------------- getLine.c --------------- */

int getLine(FILE *inputFile, char *commentChar, char *line,
	    bool_t exit_on_EOF)
{
  const char routineName[] = "getLine";

  char *linePtr, dummy[2];

  /* --- Reads (into char array line) inputFile till first line
         that does not start with commentChar and returns length of
         line, or EOF at end of file. Empty lines are also ignored.
         --                                            -------------- */

  while((linePtr = fgets(line, MAX_LINE_SIZE, inputFile)) != NULL) {
    if ((sscanf(line, "%1s", dummy) > 0) && (*line != *commentChar))
      break;
  }
  if (linePtr == NULL) {
    if (exit_on_EOF) {
      sprintf(messageStr,
	      "Reached end of input file before all data was read");
      Error(ERROR_LEVEL_2, routineName, messageStr);
    } else
      return EOF;
  } else 
    return strlen(line);
  return 0;
}
/* ------- end ---------------------------- getLine.c --------------- */

/* ------- begin -------------------------- checkNread.c ------------ */

void checkNread(int Nread, int Nrequired, const char *routineName,
		int checkPoint)
{
  /* --- Check whether Nread igreater or equal to Nrequired and issue 
         error otherwise --                           --------------- */

  if (Nread < Nrequired) {
    sprintf(messageStr, "Unable to read input file\n"
	    " At checkpoint %d. Needed at least %d, read %d item%s",
	    checkPoint, Nrequired, Nread, (Nread > 1) ? "s" : "");
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }
}
/* ------- end ---------------------------- checkNread.c ------------ */

/* ------- end ---------------------------- UpperCase.c ------------- */

void UpperCase(char *string)
{
  register int n;

  /* --- Change string content to upper case --        -------------- */

  for (n = 0;  n < (int) strlen(string);  n++)
     string[n] = toupper(string[n]);
}
/* ------- end ---------------------------- UpperCase.c ------------- */

/* ------- begin -------------------------- substring.c ------------- */

char *substring(const char *string, int N0, int Nchar)
{
  static char destination[MAX_LINE_SIZE];
  int length = strlen(string);
 
  /* --- Extract a substring of length Nchar from source string,
         starting at position N0. --                   -------------- */ 
 
  if (N0 >= length) N0 = length;
  if (N0 + Nchar >= length) Nchar = length - N0;
 
  memcpy(destination, string + N0, Nchar);
  destination[Nchar] = '\0';

  return destination;
}
/* ------- end ---------------------------- substring.c ------------- */

