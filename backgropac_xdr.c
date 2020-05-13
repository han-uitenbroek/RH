/* ------- file: -------------------------- backgropac_xdr.c --------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Tue May 12 17:52:18 2020 --

       --------------------------                      ----------RH-- */

/* --- Write separate cotributions to background opacity to file -- - */


#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "spectrum.h"
#include "background.h"
#include "error.h"
#include "inputs.h"
#include "xdr.h"


#define COMMENT_CHAR  "#"

enum opac_type {ABSORPTION, SCATTERING};


/* --- Function prototypes --                          -------------- */

bool_t write_contrib_xdr(XDR *xdrs, char *label, enum opac_type type,
			 int Nspace, double *scatt);


/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern InputData input;
extern char messageStr[];


/* ------- begin -------------------------- backgrOpac.c ------------ */

void backgrOpac(int Nlambda, double *lambda)
{
  const char routineName[] = "backgrOpac";
  register int nspect, n, k;

  char    metalID[12], format[21], contribution_out[MAX_LINE_SIZE],
          inputLine[MAX_LINE_SIZE];
  bool_t  exit_on_EOF, do_fudge;
  int     Nfudge;
  enum    opac_type type;
  double  *chi, *eta, *scatt, *thomson,
          Hmin_fudge, scatt_fudge, metal_fudge, *lambda_fudge, **fudge;
  Atom   *He;
  FILE   *fp_out, *fp_fudge;
  XDR     xdrs;

  /* --- Only evaluate for static case --              -------------- */

  atmos.moving = FALSE;
  atmos.Stokes = FALSE;
  
  if (input.solve_ne) Solve_ne(atmos.ne, TRUE);
  
  SetLTEQuantities();
  
  if (input.NonICE)
    readMolecules(MOLECULAR_CONCENTRATION_FILE);
  else
    ChemicalEquilibrium(N_MAX_CHEM_ITER, CHEM_ITER_LIMIT);

  if (strcmp(input.fudgeData, "none")) {
    do_fudge = TRUE;

    /* --- Read wavelength-dependent fudge factors to compensate for
           missing UV backround line haze --           -------------- */

    if ((fp_fudge = fopen(input.fudgeData, "r")) == NULL) {
      sprintf(messageStr, "Unable to open input file %s", input.fudgeData);
      Error(ERROR_LEVEL_2, routineName, messageStr);
    }
    sprintf(messageStr,
	    "\n-Fudging background opacities with file\n  %s\n\n",
	    input.fudgeData);
    Error(MESSAGE, routineName, messageStr);

    getLine(fp_fudge, COMMENT_CHAR, inputLine, exit_on_EOF=TRUE);
    sscanf(inputLine, "%d", &Nfudge);
    lambda_fudge = (double *) malloc(Nfudge * sizeof(double));
    fudge = matrix_double(3, Nfudge);
    for (n = 0;  n < Nfudge;  n++) {
      getLine(fp_fudge, COMMENT_CHAR, inputLine, exit_on_EOF=TRUE);
      sscanf(inputLine, "%lf %lf %lf %lf", &lambda_fudge[n],
	     &fudge[0][n], &fudge[1][n], &fudge[2][n]);
    }
    for (n = 0;  n < 3*Nfudge;  n++) fudge[0][n] += 1.0;
    fclose(fp_fudge);
  } else
    do_fudge = FALSE;

  /* --- Temporary storage for this routine --         -------------- */

  chi     = (double *) malloc(atmos.Nspace * sizeof(double));
  eta     = (double *) malloc(atmos.Nspace * sizeof(double));
  scatt   = (double *) malloc(atmos.Nspace * sizeof(double));
  thomson = (double *) malloc(atmos.Nspace * sizeof(double));

  /* --- Check whether He is present among the metals --  ----------- */

  He = (atmos.elements[1].model) ? atmos.elements[1].model : NULL;

  /* --- Go through wavelengths one-by-one --          -------------- */

  for (nspect = 0;  nspect < Nlambda;  nspect++) {

    /* --- Open output file for fractions of background opacity,
           emissivity, and scattering --               -------------- */

    sprintf(format, "bopac%%%1.1d.3f",
	    (int) log10(lambda[nspect]) + 5);
    sprintf(contribution_out, format, lambda[nspect]);
    if ((fp_out = fopen(contribution_out, "w")) == NULL) {
      sprintf(messageStr, "Unable to open output file %s",
	      contribution_out);
      Error(ERROR_LEVEL_1, routineName, messageStr);
      return;
    }
    xdrstdio_create(&xdrs, fp_out, XDR_ENCODE);

    printf("%s Processing lambda = %9.3f [nm]\n",
	   (nspect == 0) ? "\n\n " : " ", lambda[nspect]);
    xdr_double(&xdrs, &lambda[nspect]);

    /* --- Thomson scattering by free electrons --     -------------- */

    Thomson(thomson);
    if (do_fudge) {
      Linear(Nfudge, lambda_fudge, fudge[1],
	     1, &lambda[nspect], &scatt_fudge, FALSE);
      for (k = 0;  k < atmos.Nspace;  k++) thomson[k] *= scatt_fudge;
    }
    write_contrib_xdr(&xdrs, "THOMSON", type = SCATTERING,
		      atmos.Nspace, thomson);

    /* --- Negative hydrogen ion, bound-free and free-free -- ------- */

    if (do_fudge) {
      Linear(Nfudge, lambda_fudge, fudge[0],
	     1, &lambda[nspect], &Hmin_fudge, FALSE);
    }
    if (Hminus_bf(lambda[nspect], chi, eta)) {
      if (do_fudge) {
	for (k = 0;  k < atmos.Nspace;  k++) chi[k] *= Hmin_fudge;
      }
      write_contrib_xdr(&xdrs, "HMINUS_BF", type=ABSORPTION,
			atmos.Nspace, chi);
    }
    if (Hminus_ff(lambda[nspect], chi)) {
      if (do_fudge) {
	for (k = 0;  k < atmos.Nspace;  k++) chi[k] *= Hmin_fudge;
      }
      write_contrib_xdr(&xdrs, "HMINUS_FF", type=ABSORPTION,
			atmos.Nspace, chi);
    }

    /* --- Bound-free opacities from OH and CH molecules -- --------- */

    if (OH_bf_opac(lambda[nspect], chi, eta)) {
      write_contrib_xdr(&xdrs, "OH_BF", type=ABSORPTION,
			atmos.Nspace, chi);
    }
    if (CH_bf_opac(lambda[nspect], chi, eta)) {
      write_contrib_xdr(&xdrs, "CH_BF", type=ABSORPTION,
			atmos.Nspace, chi);
    }

    /* --- Neutral Hydrogen Bound-Free and Free-Free --  ------------ */

    if (Hydrogen_bf(lambda[nspect], chi, eta)) {
      write_contrib_xdr(&xdrs, "HYDROGEN_BF", type=ABSORPTION,
			atmos.Nspace, chi);
    }
    Hydrogen_ff(lambda[nspect], chi);
    write_contrib_xdr(&xdrs, "HYDROGEN_FF", type=ABSORPTION,
		      atmos.Nspace, chi);

    /* --- Rayleigh scattering by neutral hydrogen --  -------------- */

    if (Rayleigh(lambda[nspect], atmos.H, scatt)) {
      if (do_fudge) {
	for (k = 0;  k < atmos.Nspace;  k++) scatt[k] *= scatt_fudge;
      }
      write_contrib_xdr(&xdrs, "RAYLEIGH_H", type=SCATTERING,
			atmos.Nspace, scatt);
    }
    /* --- Rayleigh scattering by neutral helium --    -------------- */

    if (He && Rayleigh(lambda[nspect], He, scatt)) {
            if (do_fudge) {
	for (k = 0;  k < atmos.Nspace;  k++) scatt[k] *= scatt_fudge;
      }
      write_contrib_xdr(&xdrs, "RAYLEIGH_HE", type=SCATTERING,
			atmos.Nspace, scatt);
    }
    /* --- Rayleigh scattering and free-free absorption by
           molecular hydrogen --                       -------------- */

    if (Rayleigh_H2(lambda[nspect], scatt)) {
      if (do_fudge) {
	for (k = 0;  k < atmos.Nspace;  k++) scatt[k] *= scatt_fudge;
      }
      write_contrib_xdr(&xdrs, "RAYLEIGH_H2", type=SCATTERING,
			atmos.Nspace, scatt);
    }
    /* --- Absorption by H + H^+ (referred to as H2plus free-free) -- */

    if (H2plus_ff(lambda[nspect], chi))
      write_contrib_xdr(&xdrs, "H2PLUS_FF", type=ABSORPTION,
			atmos.Nspace, chi);

    if (H2minus_ff(lambda[nspect], chi)) {
      write_contrib_xdr(&xdrs, "H2MINUS_FF", type=ABSORPTION,
			atmos.Nspace, chi);
    }
    /* --- Bound-Free opacities due to ``metals'' --   -------------- */

    if (do_fudge) {
      Linear(Nfudge, lambda_fudge, fudge[2],
	     1, &lambda[nspect], &metal_fudge, FALSE);
    }
    for (n = 1;  n < atmos.Natom;  n++) {
      Metal_bf(lambda[nspect], 1, &atmos.atoms[n], chi, eta);
      sprintf(metalID, "METAL_BF_%2s", atmos.atoms[n].ID);
      if (do_fudge) {
	for (k = 0;  k < atmos.Nspace;  k++) chi[k] *= metal_fudge;
      }
      write_contrib_xdr(&xdrs, metalID, type=ABSORPTION,
			atmos.Nspace, chi);
    }
    xdr_destroy(&xdrs);
    fclose(fp_out);
  }

  /* --- Clean up --                                     ------------ */

  for (n = 0;  n < atmos.Natom;  n++) freeAtom(atmos.atoms + n);
  free(atmos.atoms);

  for (n = 0;  n < atmos.Nmolecule;  n++)
    freeMolecule(atmos.molecules + n);
  if (atmos.Nmolecule > 0)  free(atmos.molecules);

  /* --- Free the temporary space allocated in the ff routines -- --- */

  Hminus_ff(0.0, NULL);
  H2minus_ff(0.0, NULL);
  H2plus_ff(0.0, NULL);

  free(chi);  free(eta);  free(scatt);  free(thomson);

  if (do_fudge) {
    free(lambda_fudge);
    freeMatrix((void **) fudge);
  }
}
/* ------- end ---------------------------- backgrOpac.c ------------ */

/* ------- begin -------------------------- write_contrib_xdr.c ----- */

bool_t write_contrib_xdr(XDR *xdrs, char *label, enum opac_type type,
			 int Nspace, double *scatt)
{
  bool_t result = TRUE;

  result &= xdr_counted_string(xdrs, &label);
  result &= xdr_enum(xdrs, (enum_t *) &type);
  result &= xdr_vector(xdrs, (char *) scatt, Nspace,
		       sizeof(double), (xdrproc_t) xdr_double);

  return result;
}
/* ------- end ---------------------------- write_contrib_xdr.c ----- */
