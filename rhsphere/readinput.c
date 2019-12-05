/* ------- file: -------------------------- readinput.c -------------

       Version:       rh2.0, 1-D spherically symmetric
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Thu Jan 28 22:02:22 1999 --

       --------------------------                      ----------RH-- */

/* --- Reads input data for 1-D spherically symmetric version -- ---- */

 
#include <stdio.h>
#include <string.h>

#include "rh.h"
#include "atmos.h"
#include "spectrum.h"
#include "constant.h"
#include "statistics.h"
#include "inputs.h"
#include "error.h"


/* --- Global variables --                             -------------- */

extern struct Geometry     geometry;
extern struct Atmos        atmos;
extern struct Spectrum     spectrum;
extern struct InputData    input;
extern struct CommandLine  commandline;
extern struct ProgramStats stats;
extern char   messageStr[];


/* --- Function prototypes --                          -------------- */

void setstartValue(char *value, void *pointer);

/* ------- begin -------------------------- readInput.c ------------- */

void readInput(char *inputFileName)
{
  int   Nkeyword;
  FILE *inputFile;
  const char routineName[] = "readInput";

  struct Keyword theKeywords[] = {
    {"ATMOS_FILE", "", FALSE, KEYWORD_REQUIRED, input.atmos_input,
     setcharValue},
    {"ABUND_FILE", "", FALSE, KEYWORD_REQUIRED, input.abund_input,
     setcharValue},
    {"ATOM_FILE",  "", FALSE, KEYWORD_REQUIRED, input.atom_input,
     setcharValue},

    {"STARTING_SOLUTION", "", FALSE, KEYWORD_REQUIRED, &input.start,
     setstartValue},
    {"N_EDDINGTON", "0", FALSE, KEYWORD_OPTIONAL, &input.N_Eddington,
     setintValue}, 
    {"N_MAX_SCATTER", "5", FALSE, KEYWORD_OPTIONAL, &input.NmaxScatter,
     setintValue},

    {"I_SUM",     "0", FALSE, KEYWORD_REQUIRED, &input.iSum, setintValue},
    {"N_MAX_ITER", "", FALSE, KEYWORD_REQUIRED, &input.NmaxIter, setintValue},
    {"ITER_LIMIT", "", FALSE, KEYWORD_REQUIRED, &input.iterLimit,
     setfloatValue},
    {"NG_DELAY",  "0", FALSE, KEYWORD_OPTIONAL, &input.Ngdelay, setintValue},
    {"NG_ORDER",  "0", FALSE, KEYWORD_OPTIONAL, &input.Ngorder, setintValue},
    {"NG_PERIOD", "1", FALSE, KEYWORD_OPTIONAL, &input.Ngperiod, setintValue},
    
    {"PRD_TRESHOLD", "0.1", FALSE, KEYWORD_REQUIRED, &input.PRD_treshold,
     setfloatValue},
    {"PRD_N_MAX_ITER", "8", FALSE, KEYWORD_OPTIONAL, &input.PRD_NmaxIter,
     setintValue},
    {"PRD_NG_DELAY",  "0", FALSE, KEYWORD_DEFAULT, &input.PRD_Ngdelay,
     setintValue},
    {"PRD_NG_ORDER",  "0", FALSE, KEYWORD_DEFAULT, &input.PRD_Ngorder,
     setintValue},
    {"PRD_NG_PERIOD", "0", FALSE, KEYWORD_DEFAULT, &input.PRD_Ngperiod,
     setintValue},

    {"J_FILE",     "", FALSE, KEYWORD_REQUIRED, input.JFile, setcharValue},
    {"BACKGROUND_FILE", "", FALSE, KEYWORD_REQUIRED, input.backgroundFile,
     setcharValue},
    {"HYDROGEN_ATOM", "", FALSE, KEYWORD_REQUIRED, input.H_atom,
     setcharValue},
    {"HYDROGEN_LTE", "0", FALSE, KEYWORD_DEFAULT, &atmos.H_LTE,
     setboolValue},
    {"H2_FILE", "", FALSE, KEYWORD_REQUIRED, input.H2_molecule,
     setcharValue},
    {"KURUCZ_DATA", "none", FALSE, KEYWORD_OPTIONAL, &input.KuruczData,
     setcharValue},
    {"KURUCZ_PF_DATA", "none", FALSE, KEYWORD_OPTIONAL, &input.pfData,
     setcharValue},
    {"OPACITY_FUDGE", "none", FALSE, KEYWORD_OPTIONAL, &input.fudgeData,
     setcharValue},
    {"METALLICITY", "0.0", FALSE, KEYWORD_DEFAULT, &input.metallicity,
     setfloatValue},
    
    {"ATOM_OUTPUT",  "atom.out", FALSE, KEYWORD_DEFAULT, input.atom_output,
     setcharValue},
    {"ATMOS_OUTPUT", "atmos.out", FALSE, KEYWORD_DEFAULT, input.atmos_output,
     setcharValue},
    {"GEOMETRY_OUTPUT", "geometry.out", FALSE, KEYWORD_OPTIONAL,
     input.geometry_output, setcharValue},
    {"SPECTRUM_OUTPUT", "spectrum.out", FALSE, KEYWORD_DEFAULT,
       input.spectrum_output, setcharValue},
    {"OPACITY_OUTPUT", "none", FALSE, KEYWORD_OPTIONAL, input.opac_output,
     setcharValue},
    {"RADRATE_OUTPUT", "none", FALSE, KEYWORD_OPTIONAL, input.radrateFile,
     setcharValue},
    
    {"VMICRO_CHAR",  "",    FALSE, KEYWORD_REQUIRED, &atmos.vmicro_char,
     setfloatValue},
    {"VMACRO_TRESH", "0.1", FALSE, KEYWORD_OPTIONAL, &atmos.vmacro_tresh,
     setfloatValue},
    {"LAMBDA_REF", "500.0", FALSE, KEYWORD_DEFAULT, &atmos.lambda_ref,
     setfloatValue},
    {"VACUUM_TO_AIR", "0", FALSE, KEYWORD_OPTIONAL, &spectrum.vacuum_to_air,
     setboolValue},

    {"PRINT_CPU",       "0", FALSE, KEYWORD_OPTIONAL, &stats.printCPU,
     setboolValue}
  };
  Nkeyword = sizeof(theKeywords) / sizeof(struct Keyword);

  /* --- Open the input data file --                    ------------- */

  if ((inputFile = fopen(inputFileName, "r")) == NULL) {
    Error(UNABLE_TO_OPEN_INPUTFILE, routineName, inputFileName);
  }

  readValues(inputFile, Nkeyword, theKeywords);
  fclose(inputFile);

  /* --- Lambda_ref should be positive since it is used for the conversion
         of depthscales in routine convertScales --     ------------- */

  if (atmos.lambda_ref <= 0.0) {
    sprintf(messageStr, "Value of LAMBDA_REF should be larger than 0.0\n");
    Error(ERROR_READING_INPUTFILE, routineName, messageStr);
  } 
  /* --- Convert to MKSA units where necessary --       ------------- */

  atmos.vmicro_char  *= KM_TO_M;
  atmos.vmacro_tresh *= KM_TO_M;
}
/* ------- end ---------------------------- readInput.c ------------- */
