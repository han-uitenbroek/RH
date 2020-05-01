/* ------- file: -------------------------- rhsc2d.c ----------------

       Version:       rh2.0, 2-D Cartesian
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Fri May  1 08:25:39 2020 --

       --------------------------                      ----------RH-- */

/* --- Main routine of 2-D Cartesian radiative transfer program.
       MALI scheme formulated according to Rybicki & Hummer

  See: G. B. Rybicki and D. G. Hummer 1991, A&A 245, p. 171-181
       G. B. Rybicki and D. G. Hummer 1992, A&A 263, p. 209-215

       Formal solution is performed with the method of
       short characteristics.

  See: P. B. Kunasz and L. H. Auer 1988, JQSRT 39, p. 67

       --                                              -------------- */


#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "geometry.h"
#include "spectrum.h"
#include "background.h"
#include "statistics.h"
#include "inputs.h"
#include "xdr.h"


/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

enum Topology topology = TWO_D_PLANE;

Atmosphere atmos;
Geometry geometry;
Spectrum spectrum;
ProgramStats stats;
InputData input;
CommandLine commandline;
char messageStr[MAX_LINE_SIZE];


/* ------- begin -------------------------- rhsc2d.c ---------------- */

int main(int argc, char *argv[])
{
  bool_t analyze_output, equilibria_only;
  int    niter, nact;
  double deltaJ;

  Atom *atom;
  Molecule *molecule;

  /* --- Read input data and initialize --             -------------- */

  setOptions(argc, argv);
  getCPU(0, TIME_START, NULL);
  SetFPEtraps();

  readInput();
  spectrum.updateJ = TRUE;
 
  getCPU(1, TIME_START, NULL);
  readAtmos(&atmos, &geometry);
  if (atmos.Stokes) Bproject();
  fillMesh(&geometry);

  readAtomicModels();
  readMolecularModels();
  SortLambda();
  
  getBoundary(&atmos, &geometry);

  Background(analyze_output=TRUE, equilibria_only=FALSE);

  getProfiles();
  initSolution();
  initScatter();

  getCPU(1, TIME_POLL, "Total initialize");
 
  /* --- Solve radiative transfer for active ingredients -- --------- */

  Iterate(input.NmaxIter, input.iterLimit);

  /* --- If appropriate final solution(s) should be polarized -- ---- */
  
  adjustStokesMode(atom);

  niter = 0;
  while (niter < input.NmaxScatter) {
    deltaJ = solveSpectrum(FALSE, FALSE);
    if (!input.backgr_pol && deltaJ <= input.iterLimit) break;

    niter++;
  }

  /* --- Write output files --                     ------------------ */
 
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
/* ------- end ---------------------------- rhsc2d.c ---------------- */
