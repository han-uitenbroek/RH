/* ------- file: -------------------------- readvalues.c ------------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Mon Jun 11 14:46:13 2018 --

       --------------------------                      ----------RH-- */

/* --- Reads input parameters from file and stores them in InputData
       structure input. The input file is read line by line and
       keyword--value pairs are extracted, separated by `` = ''. Possible
       keywords are defined in the Keyword structure theKeywords. This
       structure has the following contents:

         struct Keyword {
           char keyword;               -- Strings for keyword.
           char value;                 -- String for (default) value.
           bool set;                   -- TRUE when keyword has been set
           enum keywordtype type;      -- One of three types:
                                          KEYWORD_REQUIRED, KEYWORD_DEFAULT,
                                          and KEYWORD_OPTIONAL.
           void *pointer, (*setValue)(char *value, void *pointer);
                                       -- pointer to function of type void
                                          accepting pointer to string value
                                          and pointer to actual value.
         };

       Keywords with KEYWORD_REQUIRED required always have to be set;
       a fatal error will occur otherwise. In case of a KEYWORD_DEFAULT
       the default (set in the structure initialization for theKeywords,
       in the readInput routines of each of the RH codes) is used when
       the keyword has not been set. A warning that the default has been
       overridden is issued when the keyword is set. In the case of
       KEYWORD_OPTIONAL no such warning is issued.

       The setValue pointers should point to functions that translate
       the string input value to the proper type and store the latter
       in the variable pointed to by pointer.
       --                                              -------------- */

#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "error.h"
#include "inputs.h"


#define COMMENT_CHAR "#"


/* --- Function prototypes --                          -------------- */

int pthread_setconcurrency(int new_level);


/* --- Global variables --                             -------------- */

extern InputData input;
extern CommandLine commandline;
extern char messageStr[];


/* ------- begin -------------------------- readValues.c ------------ */

void readValues(FILE *fp_keyword, int Nkeyword,	Keyword *theKeywords)
{
  const char routineName[] = "readValues";
  register int n;

  char   keyword[MAX_KEYWORD_LENGTH], value[MAX_VALUE_LENGTH],
         line[MAX_LINE_SIZE];
  bool_t recognized, exit_on_EOF;
  int    nread;

  /* --- Read input data line-by-line --                ------------- */

  while (getLine(fp_keyword, COMMENT_CHAR, line, exit_on_EOF=FALSE) != EOF) {
    if ((nread = sscanf(line, "%s = %s", keyword, value)) != 2) {
      sprintf(messageStr, "Missing input value for keyword %s", keyword);
      Error(ERROR_LEVEL_2, routineName, messageStr);
    }

    recognized = FALSE;
    for (n = 0;  n < Nkeyword;  n++) {
      if (!strcmp(keyword, theKeywords[n].keyword)) {
	recognized = TRUE;
	if (theKeywords[n].type == KEYWORD_DEFAULT) {
	  sprintf(messageStr,
		   "Overriding default value %s for keyword %s with %s",
		   theKeywords[n].value, theKeywords[n].keyword, value);
	  Error(WARNING, routineName, messageStr);
	}
        strcpy(theKeywords[n].value, value);
	theKeywords[n].set = TRUE;
        break;
      }
    }
    if (!recognized) {
      sprintf(messageStr, "Did not recognize keyword %s with value %s",
	       keyword, value);
      Error(WARNING, routineName, messageStr);
    }
  }
  /* --- Go through all possible keywords and set their values -- --- */

  for (n = 0;  n < Nkeyword;  n++) {
    if ((!theKeywords[n].set) &&
	 (theKeywords[n].type == KEYWORD_REQUIRED)) {
      sprintf(messageStr, "Missing input value for required keyword %s",
	      theKeywords[n].keyword);
      Error(ERROR_LEVEL_2, routineName, messageStr);
    }
    (*theKeywords[n].setValue)(theKeywords[n].value, theKeywords[n].pointer);
  }
}
/* ------- end ---------------------------- readValues.c ------------ */

/* ------- end ---------------------------- showValues.c ------------ */

void showValues(int Nkeyword, Keyword *theKeywords)
{
  register int n;

  char typeStr[MAX_VALUE_LENGTH], *value;
  int  length;

  fprintf(stderr, "\n  %-20s = %-30.30s %-18.18s %s\n\n",
	  "KEYWORD", "VALUE", "(KEYWORD_TYPE)", "SET");

  for (n = 0;  n < Nkeyword;  n++) {
    switch (theKeywords[n].type) {
    case KEYWORD_REQUIRED: strcpy(typeStr, "(KEYWORD_REQUIRED)");  break;
    case KEYWORD_DEFAULT:  strcpy(typeStr, "(KEYWORD_DEFAULT)");   break;
    case KEYWORD_OPTIONAL: strcpy(typeStr, "(KEYWORD_OPTIONAL)");  break;
    }
    /* --- Display end rather than begin of keyword value if string is
           too long --                                 -------------- */

    if ((length = strlen(theKeywords[n].value)) > 30)
      value = theKeywords[n].value + (length - 30);
    else
      value = theKeywords[n].value;

    fprintf(stderr, "  %-20s = %-30.30s %-18.18s  %s\n",
	    theKeywords[n].keyword, value, typeStr,
	    (theKeywords[n].set) ? "X" : "");
  }
  exit(EXIT_SUCCESS);
}
/* ------- end ---------------------------- showValues.c ------------ */

/* ------- begin -------------------------- setcharValue.c ---------- */

void setcharValue(char *value, void *pointer)
{
  strcpy((char *) pointer, value);
}
/* ------- end ---------------------------- setcharValue.c ---------- */

/* ------- begin -------------------------- setintValue.c ----------- */

void setintValue(char *value, void *pointer)
{
  int intvalue = atoi(value);

  memcpy(pointer, &intvalue, sizeof(int));
}
/* ------- end ---------------------------- setintValue.c ----------- */

/* ------- begin -------------------------- setdoubleValue.c -------- */

void setdoubleValue(char *value, void *pointer)
{
  double doublevalue = atof(value);

  memcpy(pointer, &doublevalue, sizeof(double));
}
/* ------- end ---------------------------- setdoubleValue.c -------- */

/* ------- begin -------------------------- setboolValue.c ---------- */

void setboolValue(char *value, void *pointer)
{
  bool_t boolvalue = (strstr(value, "TRUE")) ? TRUE : FALSE;

  memcpy(pointer, &boolvalue, sizeof(bool_t));
}
/* ------- end ---------------------------- setboolValue.c ---------- */

/* ------- begin -------------------------- setAngleSet.c ----------- */

void setAngleSet(char *value, void *pointer)
{
  int  Nread, Ninclination, Nazimuth;
  AngleSet angleSet;

  initAngleSet(&angleSet);

  if (!strcmp(value, "SET_VERTICAL")) angleSet.set = SET_VERTICAL;
  else if  (!strcmp(value, "SET_A2")) angleSet.set = SET_A2;
  else if  (!strcmp(value, "SET_A4")) angleSet.set = SET_A4;
  else if  (!strcmp(value, "SET_A6")) angleSet.set = SET_A6;
  else if  (!strcmp(value, "SET_A8")) angleSet.set = SET_A8;
  else if  (!strcmp(value, "SET_B4")) angleSet.set = SET_B4;
  else if  (!strcmp(value, "SET_B6")) angleSet.set = SET_B6;
  else if  (!strcmp(value, "SET_B8")) angleSet.set = SET_B8;
  else if  (!strcmp(value, "NO_SET")) angleSet.set = NO_SET;

  else if (strstr(value, "SET_GL_")) {
    if (sscanf(value, "SET_GL_%dX%d", &Ninclination, &Nazimuth) != 2) {
      sprintf(messageStr,
	      "\n  Invalid Gauss-Legandre format for keyword ANGLE_SET: %s",
	      value);
      Error(ERROR_LEVEL_2, "setAngleSet", messageStr);
    }
    angleSet.set = SET_GL;
    angleSet.Ninclination = Ninclination;
    angleSet.Nazimuth = Nazimuth;

    sprintf(messageStr,
	    "\n  Found Gauss-Legendre angleset with "
	    "Ninclination = %d and Nazimuth = %d", 
	    Ninclination, Nazimuth);
    Error(MESSAGE, "setAngleSet", messageStr);

    if (Ninclination <= 0  || Ninclination > NMAXINCLINATION) {
      sprintf(messageStr,
	      "\n  Invalid value for Ninclination in angleset: %d",
	      Ninclination);
      Error(ERROR_LEVEL_2, "setAngleSet", messageStr);
    }
    if (Nazimuth <= 0  || Nazimuth > NMAXAZIMUTH) {
      sprintf(messageStr,
	      "\n  Invalid value for Nazimuth in angleset: %d",
	      Nazimuth);
      Error(ERROR_LEVEL_2, "setAngleSet", messageStr);
    }
  } else {
    sprintf(messageStr,
	    "\n  Invalid value for keyword ANGLE_SET: %s", value);
    Error(ERROR_LEVEL_2, "setAngleSet", messageStr);
  }
  memcpy(pointer, &angleSet, sizeof(AngleSet));
}
/* ------- begin -------------------------- initAngleSet.c ---------- */

void initAngleSet(AngleSet *angleSet)
{
  angleSet->set = NO_SET;
  angleSet->Ninclination = 0;
  angleSet->Nazimuth = 0;
}
/* ------- end ---------------------------- initAngleSet.c ---------- */

/* ------- begin -------------------------- setStokesMode.c --------- */

void setStokesMode(char *value, void *pointer)
{
  const char routineName[] = "setStokesMode";

  enum StokesMode StokesMode;

  if (!strcmp(value, "NO_STOKES"))
    StokesMode = NO_STOKES;
  else if (!strcmp(value, "FIELD_FREE"))
    StokesMode = FIELD_FREE;
  else if (!strcmp(value, "POLARIZATION_FREE"))
    StokesMode = POLARIZATION_FREE;
  else if (!strcmp(value, "FULL_STOKES"))
    StokesMode = FULL_STOKES;
  else {
    sprintf(messageStr,
	     "Invalid value for keyword STOKES_MODE: %s", value);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }
  memcpy(pointer, &StokesMode, sizeof(enum_t));
}
/* ------- end ---------------------------- setStokesMode.c --------- */

/* ------- begin -------------------------- setThreadValue.c -------- */

#define N_THREAD_LIMIT 32

void setThreadValue(char *value, void *pointer)
{
  const char routineName[] = "setThreadValue";

  int Nthreads = atoi(value), return_value;

  if (Nthreads > N_THREAD_LIMIT) {
    sprintf(messageStr,
	    "Value of keyword N_THREADS (%d) larger than allowed limit (%d)",
	    Nthreads, N_THREAD_LIMIT);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  } else if (Nthreads < 1) {
    Nthreads = 1;
  }

  if (Nthreads > 1) {

    /* --- For now we try to use the default thread attributes,
           most notably:

          - PTHREAD_SCOPE_PROCESS   for scope
          - PTHREAD_CREATE_JOINABLE for detach state
          - PTHREAD_INHERIT_SCHED   for scheduling inheritance
          - SCHED_OTHER             for scheduling policy
	  ---                                          -------------- */

    pthread_attr_init(&input.thread_attr);

    if (pthread_attr_setscope(&input.thread_attr, PTHREAD_SCOPE_PROCESS))
      Error(WARNING, routineName, "Non-default thread scope");
    if (pthread_attr_setdetachstate(&input.thread_attr,
				    PTHREAD_CREATE_JOINABLE))
      Error(WARNING, routineName, "Non-default thread detach state");
    if (pthread_attr_setinheritsched(&input.thread_attr,
				     PTHREAD_INHERIT_SCHED))
      Error(WARNING, routineName, "Non-default thread scheduling inheritance");
    if (pthread_attr_setschedpolicy(&input.thread_attr, SCHED_OTHER))
      Error(WARNING, routineName, "Non-default thread scheduling policy");

    /* --- Suggest thread concurrency to the operating system -- ---- */
    
    if ((return_value = pthread_setconcurrency(Nthreads)))
      Error(ERROR_LEVEL_2, routineName,
	    "Failed to set concurrency level for threads.");
    else {
      sprintf(messageStr, "Setting thread concurrency to %d", Nthreads);
      Error(WARNING, routineName, messageStr);
    }
  }

  memcpy(pointer, &Nthreads, sizeof(int));
}
/* ------- end ---------------------------- setThreadValue.c -------- */

/* ------- begin -------------------------- setInterpolate_3D.c ----- */

void setInterpolate_3D(char *value, void *pointer)
{
  const char routineName[] = "setInterpolate_3D";

  enum order_3D order;

  if (!strcmp(value, "LINEAR_3D"))
    order = LINEAR_3D;
  else if (!strcmp(value, "BICUBIC_3D"))
    order = BICUBIC_3D;
  else {
    sprintf(messageStr,
	    "\n  Invalid value for keyword INTERPOLATE_3D: %s", value);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }
  memcpy(pointer, &order, sizeof(enum order_3D));
}
/* ------- end ---------------------------- setInterpolate_3D.c ----- */


/* ------- begin -------------------------- set_S_Interpolation.c --- */

void set_S_Interpolation(char *value, void *pointer)
{
  const char routineName[] = "setInterpolate_3D";

  enum S_interpol interpolation;

  if (!strcmp(value, "S_PARABOLIC"))
    interpolation = S_PARABOLIC;
  else if (!strcmp(value, "S_LINEAR"))
    interpolation = S_LINEAR;
  else if (!strcmp(value, "S_BEZIER3"))
    interpolation = S_BEZIER3;
  else {
    sprintf(messageStr,
	    "\n  Invalid value for keyword S_INTERPOLATION: %s", value);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }
  memcpy(pointer, &interpolation, sizeof(enum S_interpol));
}
/* ------- end ---------------------------- set_S_Interpolation.c --- */


/* ------- begin ------------------- set_S_interpolation_stokes.c --- */

void set_S_interpolation_stokes(char *value, void *pointer)
{
  const char routineName[] = "set_S_interpolation_stokes";

  enum S_interpol_stokes interpolation;

  if (!strcmp(value, "DELO_BEZIER3"))
    interpolation = DELO_BEZIER3;
  else if (!strcmp(value, "DELO_PARABOLIC"))
    interpolation = DELO_PARABOLIC;
  else {
    sprintf(messageStr,
	     "Invalid value for keyword S_INTERPOLATION_STOKES: %s", value);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }
  memcpy(pointer, &interpolation, sizeof(enum S_interpol_stokes));
}
/* ------- end --------------------- set_S_interpolation_stokes.c --- */


/* ------- begin -------------------------- setstartValue.c --------- */

void setstartValue(char *value, void *pointer)
{
  const char routineName[] = "setStartValue";

  enum solution startvalue;

  if (!strcmp(value, "NEW_J"))
    startvalue = NEW_J;
  else if (!strcmp(value, "OLD_J"))
    startvalue = OLD_J;
  else {
    sprintf(messageStr,
             "Invalid value for keyword STARTING_J: %s", value);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }

  memcpy(pointer, &startvalue, sizeof(enum_t));
}
/* ------- end ---------------------------- setstartValue.c --------- */

/* ------- begin -------------------------- setnesolution.c --------- */

void setnesolution(char *value, void *pointer)
{
  const char routineName[] = "setnesolution";

  enum ne_solution nesolution;

  if (!strcmp(value, "NONE"))
    nesolution = NONE;
  else if (!strcmp(value, "ONCE"))
    nesolution = ONCE;
  else if (!strcmp(value, "ITERATION")) {
    nesolution = ITERATION;
  } else {
    sprintf(messageStr,
             "Invalid value for keyword SOLVE_NE: %s", value);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }

  memcpy(pointer, &nesolution, sizeof(enum_t));
}
/* ------- end ---------------------------- setnesolution.c --------- */
