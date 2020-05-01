/* ------- file: -------------------------- rhsc3d.c ----------------

       Version:       rh2.0, 3-D Cartesian, short characteristics
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Fri May  1 10:07:12 2020 --

       --------------------------                      ----------RH-- */

/* --- Main routine of 3-D Cartesian radiative transfer program.
       MALI scheme formulated according to Rybicki & Hummer.

  See: G. B. Rybicki and D. G. Hummer 1991, A&A 245, p. 171-181
       G. B. Rybicki and D. G. Hummer 1992, A&A 263, p. 209-215

       Formal solution is performed with the short characteristics
       method.

  See: P. B. Kunasz and L. H. Auer 1988, JQSRT 39, p. 67

       Interpolation in the horizontal planes is done with cubic
       convolution.

  See: R.G. Keys, 1981, in IEEE Trans. Acoustics, Speech,
       and Signal Processing, Vol. 29, pp. 1153-1160.

       The horizontal grid should be equidistant and periodic in both
       x- and y directions, although not necessarily with the same
       period or separation in x and y.

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


enum Topology topology = THREE_D_PLANE;

Atmosphere atmos;
Geometry geometry;
Spectrum spectrum;
ProgramStats stats;
InputData input;
CommandLine commandline;
char messageStr[MAX_LINE_SIZE];


/* ------- begin -------------------------- main.c ------------------ */

int main(int argc, char *argv[])
{
  bool_t analyze_output, equilibria_only;
  int    niter, nact;
  double deltaJ;
  
  Atom *atom;
  Molecule *molecule;

  setOptions(argc, argv);
  getCPU(0, TIME_START, NULL);
  SetFPEtraps();

  /* --- Read input data and initialize --             -------------- */

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

  adjustStokesMode();
  
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
/* ------- end ---------------------------- main.c ------------------ */
