/* ------- file: -------------------------- rhf1d.c -----------------

       Version:       rh2.0, 1-D plane-parallel
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Thu Feb 24 16:40:14 2011 --

       --------------------------                      ----------RH-- */

/* --- Main routine of 1D plane-parallel radiative transfer program.
       MALI scheme formulated according to Rybicki & Hummer

  See: G. B. Rybicki and D. G. Hummer 1991, A&A 245, p. 171-181
       G. B. Rybicki and D. G. Hummer 1992, A&A 263, p. 209-215

       Formal solution is performed with Feautrier difference scheme
       in static atmospheres, and with piecewise quadratic integration
       in moving atmospheres.

       --                                              -------------- */

#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "geometry.h"
#include "spectrum.h"
#include "background.h"
#include "statistics.h"
#include "error.h"
#include "inputs.h"
#include "xdr.h"


/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

enum Topology topology = ONE_D_PLANE;

Atmosphere atmos;
Geometry geometry;
Spectrum spectrum;
ProgramStats stats;
InputData input;
CommandLine commandline;
char messageStr[MAX_MESSAGE_LENGTH];


/* ------- begin -------------------------- rhf1d.c ----------------- */

int main(int argc, char *argv[])
{
  bool_t write_analyze_output, equilibria_only;
  int    niter, nact;

  Atom *atom;
  Molecule *molecule;

  /* --- Read input data and initialize --             -------------- */

  setOptions(argc, argv);
  getCPU(0, TIME_START, NULL);
  SetFPEtraps();

  readInput();
  spectrum.updateJ = TRUE;

  getCPU(1, TIME_START, NULL);
  MULTIatmos(&atmos, &geometry);
  if (atmos.Stokes) Bproject();

  readAtomicModels();
  readMolecularModels();
  SortLambda();
  
  getBoundary(&geometry);
  
  Background(write_analyze_output=TRUE, equilibria_only=FALSE);
  convertScales(&atmos, &geometry);

  getProfiles();
  initSolution();
  initScatter();

  getCPU(1, TIME_POLL, "Total Initialize");

  /* --- Solve radiative transfer for active ingredients -- --------- */

  Iterate(input.NmaxIter, input.iterLimit);

  adjustStokesMode();
  niter = 0;
  while (niter < input.NmaxScatter) {
    if (solveSpectrum(FALSE, FALSE) <= input.iterLimit) break;
    niter++;
  }
  /* --- Write output files --                         -------------- */

  if (atmos.hydrostatic) {
    geometry.scale = COLUMN_MASS;
    convertScales(&atmos, &geometry);
  }
  getCPU(1, TIME_START, NULL);

  writeInput();
  writeAtmos(&atmos);
  writeGeometry(&geometry);
  writeSpectrum(&spectrum);
  writeFlux(FLUX_DOT_OUT);

  for (nact = 0;  nact < atmos.Nactiveatom;  nact++) {
    atom = atmos.activeatoms[nact];

    writeAtom(atom);
    writePopulations(atom);
    writeRadRate(atom);
    writeCollisionRate(atom);
    writeDamping(atom);
  } 
  for (nact = 0;  nact < atmos.Nactivemol;  nact++) {
    molecule = atmos.activemols[nact];
    writeMolPops(molecule);
  }

  writeOpacity();

  getCPU(1, TIME_POLL, "Write output");
  printTotalCPU();
}
/* ------- end ---------------------------- rhf1d.c ----------------- */
