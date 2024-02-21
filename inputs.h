/* ------- file: -------------------------- inputs.h ----------------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Tue Feb 20 09:34:35 2024 --

       --------------------------                      ----------RH-- */

#ifndef __INPUTS_H__
#define __INPUTS_H__


/* --- Include file for command line options and keyword input
       parameters. --                                  -------------- */


#define MAX_KEYWORD_LENGTH  32
#define MAX_VALUE_LENGTH    160

#define CHECK_OPTION(option, name, length) \
        (strncmp(option, name, MAX(((int) length), strlen(option))) == 0)


enum keywordtype  {KEYWORD_REQUIRED, KEYWORD_DEFAULT, KEYWORD_OPTIONAL};

enum S_interpol          {S_LINEAR, S_PARABOLIC, S_BEZIER3};
enum S_interpol_stokes   {DELO_PARABOLIC, DELO_BEZIER3};

enum order_3D     {LINEAR_3D, BICUBIC_3D};
enum ne_solution  {NONE, ONCE, ITERATION};


typedef struct {
  char   name[MAX_KEYWORD_LENGTH];
  int    minlength;
  bool_t value_required;
  char   value[MAX_VALUE_LENGTH];
  void  *pointer, (*setValue)(char *value, void *pointer);
  char   message[MAX_VALUE_LENGTH];
} Option;

typedef struct {
  char   keyword_input[MAX_VALUE_LENGTH], *wavetable, *molecule;
  bool_t quiet, showkeywords;
  FILE  *logfile;
} CommandLine;

typedef struct {
  char   keyword[MAX_KEYWORD_LENGTH], value[MAX_VALUE_LENGTH];
  bool_t set;
  enum   keywordtype type;
  void  *pointer, (*setValue)(char *value, void *pointer);
} Keyword;

typedef struct {
  char   atmos_input[MAX_VALUE_LENGTH],
         abund_input[MAX_VALUE_LENGTH],
         wavetable_input[MAX_VALUE_LENGTH],
         atoms_input[MAX_VALUE_LENGTH],
         molecules_input[MAX_VALUE_LENGTH],
         Stokes_input[MAX_VALUE_LENGTH],
         KuruczData[MAX_VALUE_LENGTH],
         pfData[MAX_VALUE_LENGTH],
         fudgeData[MAX_VALUE_LENGTH],
         atmos_output[MAX_VALUE_LENGTH],
         spectrum_output[MAX_VALUE_LENGTH],
         geometry_output[MAX_VALUE_LENGTH],
         opac_output[MAX_VALUE_LENGTH],
         JFile[MAX_VALUE_LENGTH],
         background_File[MAX_VALUE_LENGTH],
         background_ray_File[MAX_VALUE_LENGTH],
         H_atom[MAX_VALUE_LENGTH],
         H2_molecule[MAX_VALUE_LENGTH],
         radrateFile[MAX_VALUE_LENGTH],
         collrateFile[MAX_VALUE_LENGTH],
         dampingFile[MAX_VALUE_LENGTH],
         coolingFile[MAX_VALUE_LENGTH],
         Itop[MAX_VALUE_LENGTH];
  bool_t magneto_optical, XRD, Eddington,
         backgr_pol, limit_memory, allow_passive_bb, NonICE,
         rlkscatter, xdr_endian, old_background, accelerate_mols,
         RLK_explicit, prdh_limit_mem;
  enum   solution startJ;
  
  enum   StokesMode StokesMode;
  enum   S_interpol S_interpolation;
  enum   S_interpol_stokes S_interpolation_stokes;
  enum   PRDangle PRD_angle_dep;
  enum   order_3D interpolate_3D;
  enum   ne_solution solve_ne;
  
  int    isum, Ngdelay, Ngorder, Ngperiod, NmaxIter,
         PRD_NmaxIter, PRD_Ngdelay, PRD_Ngorder, PRD_Ngperiod,
         NmaxScatter, Nthreads, CR_Nstep;
  double iterLimit, PRDiterLimit, metallicity, CR_factor;

  pthread_attr_t thread_attr;
} InputData;


/* --- Associated function prototypes --               -------------- */

int   getLine(FILE *inputFile, char *commentChar, char *line,
	      bool_t exit_on_EOF);
void  parse(int argc, char *argv[], int Noption, Option *theOptions);
void  readInput();
void  readValues(FILE *fp_keyword, int Nkeyword, Keyword *theKeywords);

void  setAngleSet(char *value, void *pointer);
void  setcharValue(char *value, void *pointer);
void  setboolValue(char *value, void *pointer);
void  setintValue(char *value, void *pointer);
void  setdoubleValue(char *value, void *pointer);
void  setstartValue(char *value, void *pointer);
void  setnesolution(char *value, void *pointer);
void  setStokesMode(char *value, void *pointer);
void  setPRDangle(char *value, void *pointer);
void  setThreadValue(char *value, void *pointer);

void  setInterpolate_3D(char *value, void *pointer);
void  set_S_Interpolation(char *value, void *pointer);
void  set_S_interpolation_stokes(char *value, void *pointer);

void  showValues(int Nkeyword, Keyword *theKeywords);
char *substring(const char *string, int N0, int Nchar);
void  UpperCase(char *string);

void  writeInput();

#endif /* !__INPUTS_H__ */

/* ------- end ---------------------------- inputs.h ---------------- */
