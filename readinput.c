/* ------- file: -------------------------- readinput.c -------------

       Version:       rh2.0, 1-D plane-parallel
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Wed Oct  9 16:26:55 2019 --

       --------------------------                      ----------RH-- */

/* --- Reads input data for and defines keywords for 1-D
       plane-parallel version --                       -------------- */

 
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "spectrum.h"
#include "constant.h"
#include "statistics.h"
#include "inputs.h"
#include "error.h"


/* --- Global variables --                             -------------- */

extern enum Topology topology;

extern Atmosphere atmos;
extern Spectrum spectrum;
extern InputData input;
extern CommandLine commandline;
extern ProgramStats stats;
extern char messageStr[];


/* --- Function prototypes --                          -------------- */


/* ------- begin -------------------------- readInput.c ------------- */

void readInput()
{
  const char routineName[] = "readInput";
  static char atom_input[MAX_VALUE_LENGTH], molecule_input[MAX_VALUE_LENGTH];

  int   Nkeyword;
  FILE *fp_keyword;

  Keyword theKeywords[] = {
    {"ATMOS_FILE", "", FALSE, KEYWORD_REQUIRED, input.atmos_input,
     setcharValue},
    {"ABUND_FILE", "", FALSE, KEYWORD_REQUIRED, input.abund_input,
     setcharValue},

    {"NRAYS",     "0", FALSE, KEYWORD_OPTIONAL, &atmos.Nrays, setintValue},
    {"ANGLE_SET", "NO_SET", FALSE, KEYWORD_OPTIONAL, &atmos.angleSet,
     setAngleSet},

    {"EDDINGTON", "FALSE", FALSE, KEYWORD_OPTIONAL, &input.Eddington,
     setboolValue},
    {"ATMOS_ITOP", "none", FALSE, KEYWORD_OPTIONAL, input.Itop, setcharValue},

    {"WAVETABLE", "none", FALSE, KEYWORD_OPTIONAL, input.wavetable_input,
     setcharValue},
    {"ATOMS_FILE",  "", FALSE, KEYWORD_REQUIRED, input.atoms_input,
     setcharValue},
    {"MOLECULES_FILE",  "", FALSE, KEYWORD_REQUIRED, input.molecules_input,
     setcharValue},
    {"NON_ICE", "FALSE", FALSE, KEYWORD_OPTIONAL, &input.NonICE,
     setboolValue},

    {"N_MAX_SCATTER", "5", FALSE, KEYWORD_OPTIONAL, &input.NmaxScatter,
     setintValue},

    {"I_SUM",     "0", FALSE, KEYWORD_REQUIRED, &input.isum, setintValue},
    {"N_MAX_ITER", "", FALSE, KEYWORD_REQUIRED, &input.NmaxIter,
     setintValue},
    {"ITER_LIMIT", "", FALSE, KEYWORD_REQUIRED, &input.iterLimit,
     setdoubleValue},
    {"NG_DELAY",  "0", FALSE, KEYWORD_OPTIONAL, &input.Ngdelay, setintValue},
    {"NG_ORDER",  "0", FALSE, KEYWORD_OPTIONAL, &input.Ngorder, setintValue},
    {"NG_PERIOD", "1", FALSE, KEYWORD_OPTIONAL, &input.Ngperiod,
     setintValue},
    {"NG_MOLECULES", "FALSE", FALSE, KEYWORD_DEFAULT, &input.accelerate_mols,
     setboolValue},
    {"PRD_N_MAX_ITER", "3", FALSE, KEYWORD_OPTIONAL, &input.PRD_NmaxIter,
     setintValue},
    {"PRD_ITER_LIMIT", "1.0E-2", FALSE, KEYWORD_OPTIONAL, &input.PRDiterLimit,
     setdoubleValue},
    {"PRD_NG_DELAY",  "0", FALSE, KEYWORD_OPTIONAL, &input.PRD_Ngdelay,
     setintValue},
    {"PRD_NG_ORDER",  "0", FALSE, KEYWORD_OPTIONAL, &input.PRD_Ngorder,
     setintValue},
    {"PRD_NG_PERIOD", "0", FALSE, KEYWORD_OPTIONAL, &input.PRD_Ngperiod,
     setintValue},
    {"PRD_ANGLE_DEP", "FALSE", FALSE, KEYWORD_DEFAULT, &input.PRD_angle_dep,
     setboolValue},
    {"XRD", "FALSE", FALSE, KEYWORD_DEFAULT, &input.XRD, setboolValue}, 

    {"J_FILE",     "", FALSE, KEYWORD_REQUIRED, input.JFile, setcharValue},
    {"BACKGROUND_FILE", "", FALSE, KEYWORD_REQUIRED, input.background_File,
     setcharValue},
    {"BACKGROUND_RAY_FILE", "background.ray", FALSE, 
     KEYWORD_OPTIONAL, input.background_ray_File,
     setcharValue},
    {"OLD_BACKGROUND", "FALSE", FALSE, KEYWORD_OPTIONAL,
     &input.old_background, setboolValue},
    {"STARTING_J", "", FALSE, KEYWORD_REQUIRED, &input.startJ,
     setstartValue},
    {"HYDROGEN_LTE", "FALSE", FALSE, KEYWORD_DEFAULT, &atmos.H_LTE,
     setboolValue},
    {"HYDROSTATIC", "FALSE", FALSE, KEYWORD_DEFAULT, &atmos.hydrostatic,
     setboolValue},
    {"KURUCZ_DATA", "none", FALSE, KEYWORD_OPTIONAL, &input.KuruczData,
     setcharValue},
    {"RLK_SCATTER", "FALSE", FALSE, KEYWORD_DEFAULT, &input.rlkscatter,
     setboolValue},
    {"KURUCZ_PF_DATA", "../../Atoms/pf_Kurucz.input", FALSE,
     KEYWORD_REQUIRED, &input.pfData, setcharValue},
    {"SOLVE_NE", "NONE", FALSE, KEYWORD_DEFAULT, &input.solve_ne,
     setnesolution},
    {"OPACITY_FUDGE", "none", FALSE, KEYWORD_OPTIONAL, &input.fudgeData,
     setcharValue},
    {"METALLICITY", "0.0", FALSE, KEYWORD_DEFAULT, &input.metallicity,
     setdoubleValue},
    
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
    {"COLLRATE_OUTPUT", "none", FALSE, KEYWORD_OPTIONAL, input.collrateFile,
     setcharValue},
    {"DAMPING_OUTPUT", "none", FALSE, KEYWORD_OPTIONAL, input.dampingFile,
     setcharValue},
    {"COOLING_OUTPUT", "none", FALSE, KEYWORD_OPTIONAL, input.coolingFile,
     setcharValue},

    {"VMICRO_CHAR", "",     FALSE, KEYWORD_REQUIRED, &atmos.vmicro_char,
     setdoubleValue},
    {"VMACRO_TRESH", "0.1", FALSE, KEYWORD_OPTIONAL, &atmos.vmacro_tresh,
     setdoubleValue},
    {"LAMBDA_REF",   "500.0", FALSE, KEYWORD_DEFAULT, &atmos.lambda_ref,
     setdoubleValue},
    {"VACUUM_TO_AIR", "0", FALSE, KEYWORD_OPTIONAL, &spectrum.vacuum_to_air,
     setboolValue},

    /* --- Magnetic field related inputs go here --     ------------- */

    {"STOKES_INPUT", "none", FALSE, KEYWORD_OPTIONAL, input.Stokes_input,
     setcharValue},
    {"B_STRENGTH_CHAR", "0.0", FALSE, KEYWORD_DEFAULT, &atmos.B_char,
     setdoubleValue},
    {"STOKES_MODE", "NO_STOKES", FALSE, KEYWORD_OPTIONAL,
     &input.StokesMode, setStokesMode},
    {"MAGNETO_OPTICAL", "TRUE", FALSE, KEYWORD_DEFAULT,
     &input.magneto_optical, setboolValue},
    {"BACKGROUND_POLARIZATION", "FALSE", FALSE, KEYWORD_DEFAULT,
     &input.backgr_pol, setboolValue},
    {"XDR_ENDIAN", "TRUE", FALSE, KEYWORD_OPTIONAL,
     &input.xdr_endian, setboolValue},

    {"S_INTERPOLATION", "S_BEZIER3", FALSE, KEYWORD_DEFAULT,
     &input.S_interpolation, set_S_Interpolation},
    {"S_INTERPOLATION_STOKES", "DELO_PARABOLIC", FALSE, KEYWORD_DEFAULT,
     &input.S_interpolation_stokes, set_S_interpolation_stokes},

    {"INTERPOLATE_3D", "BICUBIC_3D", FALSE, KEYWORD_DEFAULT,
     &input.interpolate_3D, setInterpolate_3D},
    
    {"PRINT_CPU", "0", FALSE, KEYWORD_OPTIONAL, &stats.printCPU,
     setboolValue},
    {"N_THREADS", "0", FALSE, KEYWORD_OPTIONAL, &input.Nthreads,
     setThreadValue},

    {"LIMIT_MEMORY", "FALSE", FALSE, KEYWORD_DEFAULT, &input.limit_memory,
     setboolValue},
    {"ALLOW_PASSIVE_BB", "TRUE", FALSE, KEYWORD_DEFAULT,
     &input.allow_passive_bb, setboolValue}
  };
  Nkeyword = sizeof(theKeywords) / sizeof(Keyword);

  /* --- Open the input data file --                    ------------- */

  if ((fp_keyword = fopen(commandline.keyword_input, "r")) == NULL) {
    sprintf(messageStr, "Unable to open inputfile %s",
	    commandline.keyword_input);
    Error(ERROR_LEVEL_2, routineName, messageStr);
  }

  readValues(fp_keyword, Nkeyword, theKeywords);
  fclose(fp_keyword);

  /* --- Perform some sanity checks --                 -------------- */

  switch (topology) {
  case ONE_D_PLANE:
    if (atmos.Nrays == 0) {
      Error(ERROR_LEVEL_2, routineName,
	    "Must set keyword NRAYS in 1-D plane parallel geometry");
    }
    if (atmos.angleSet.set != NO_SET) {
      Error(WARNING, routineName,
	    "Ignoring value of keyword ANGLE_SET in 1-D plane geometry");
    }
    break;

  case SPHERICAL_SYMMETRIC:
    if (atmos.Nrays > 0 || atmos.angleSet.set != NO_SET) {
      Error(WARNING, routineName,
	    "Ignoring value of keywords ANGLE_SET and NRAYS in "
            "spherical geometry");
    }
    break;
 
  case TWO_D_PLANE:
    if (input.Eddington && atmos.angleSet.set != NO_SET) {
      Error(WARNING, routineName,
	    "Ignoring value of keywords ANGLE_SET > NO_VERTICAL when\n "
	    " EDDINGTON is set to TRUE\n"
	    " Using SET_VERTICAL with muz = 1/sqrt(3)");
      atmos.angleSet.set = SET_VERTICAL;
    }
    break;

  case THREE_D_PLANE:
    if (atmos.angleSet.set == NO_SET)
      Error(ERROR_LEVEL_2, routineName,
	    "Must set keyword ANGLE_SET in multi-D plane geometry");
    if (atmos.Nrays > 0)
      Error(WARNING, routineName,
	    "Ignoring value of keyword NRAYS in multi-D plane geometry");
    break;
  }

  /* --- Lambda_ref should be positive since it is used for the conversion
         of depth scales in routine convertScales --    ------------- */

  switch (topology) {
  case ONE_D_PLANE:
  case SPHERICAL_SYMMETRIC:
    if (atmos.lambda_ref <= 0.0)
      Error(ERROR_LEVEL_2, routineName,
	    "Value of LAMBDA_REF should be larger than 0.0");
    break;
  default:
    if (atmos.lambda_ref < 0.0)
      Error(ERROR_LEVEL_2, routineName,
	    "Value of LAMBDA_REF should be larger than or equal 0.0");
    break;
  }
  /* --- Stokes for the moment only in 1D plane --     -------------- */
 
  if (strcmp(input.Stokes_input, "none")) {
    switch (topology) {
    case ONE_D_PLANE:
    case TWO_D_PLANE:
      if (atmos.B_char == 0.0) {
	Error(WARNING, routineName,
	      "Parameter atmos.B_char not set or set to zero\n"
	      " Wavelength grids of line profiles do not take account "
	      " of Zeeman splitting");
      }
      if (input.StokesMode == NO_STOKES) {
        sprintf(messageStr, "%s",
	      "Keyword STOKES_MODE == NO_STOKES.\n"
	      " Set to FIELD_FREE, POLARIZATION_FREE, or FULL_STOKES\n"
	      " when doing polarization calculations");
	Error(ERROR_LEVEL_1, routineName, messageStr);
      }
      break;
    case THREE_D_PLANE:
      if (atmos.B_char == 0.0) {
	Error(WARNING, routineName,
	      "Parameter atmos.B_char not set or set to zero\n"
	      " Wavelength grids of line profiles do not take account "
	      " of Zeeman splitting");
      }
      if (input.StokesMode == NO_STOKES) {
        sprintf(messageStr, "%s",
	      "Keyword STOKES_MODE == NO_STOKES.\n"
	      " Set to FIELD_FREE, POLARIZATION_FREE, or FULL_STOKES\n"
	      " when doing polarization calculations");
	Error(ERROR_LEVEL_1, routineName, messageStr);
      }
      break;
    case SPHERICAL_SYMMETRIC:
      Error(ERROR_LEVEL_2, routineName,
	    "Cannot accomodate magnetic fields in this topology");
    }
  } else {
    if (atmos.B_char != 0.0) {
      Error(WARNING, routineName,
	    "Ignoring value of keyword B_STRENGTH_CHAR when no "
	    "magnetic field is read");
    }
    if (input.StokesMode > NO_STOKES) {
      Error(WARNING, routineName,
	    "Ignoring value of keyword STOKES_MODE when no "
	    "magnetic field is read");
    }
  }
  /* --- Hydrostatic equilibrium only in 1-D plane parallel -- ------ */

  switch (topology) {
  case ONE_D_PLANE:
    break;
  default:
    if (atmos.hydrostatic)
      Error(ERROR_LEVEL_2, routineName,
	    "Can only perform hydrostatic equilibrium calculation"
            " in 1-D Cartesian geometry");
    break;
  }
  
  /* --- If called with -showkeywords commandline option -- --------- */

  if (commandline.showkeywords) showValues(Nkeyword, theKeywords);

  /* --- Convert to MKSA units where necessary --       ------------- */

  atmos.vmicro_char  *= KM_TO_M;
  atmos.vmacro_tresh *= KM_TO_M;
}
/* ------- end ---------------------------- readInput.c ------------- */
