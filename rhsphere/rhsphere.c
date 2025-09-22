/* ------- file: -------------------------- rhspere.c ---------------

       Version:       rh2.0, 1-D spherically symmetric
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Mon Sep 22 16:41:03 2025 --

       --------------------------                      ----------RH-- */

/* --- Main routine of 1-D spherical-symmetric radiative transfer program.
       MALI scheme formulated according to Rybicki & Hummer
 
  See: G. B. Rybicki and D. G. Hummer 1991, A&A 245, p. 171-181
       G. B. Rybicki and D. G. Hummer 1992, A&A 263, p. 209-215
 
       Formal solution is performed with Feautrier difference scheme
 
       --                                              -------------- */
  
#include <stdlib.h>
#include <math.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "geometry.h"
#include "spectrum.h"
#include "background.h"
#include "inputs.h"
#include "statistics.h"
#include "xdr.h"


/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

enum Topology topology = SPHERICAL_SYMMETRIC;

Atmosphere atmos;
Geometry geometry;
Spectrum spectrum;
InputData input;
ProgramStats stats;
CommandLine commandline;
char messageStr[MAX_LINE_SIZE];


/* ------- begin -------------------------- rhsphere.c -------------- */

int main(int argc, char *argv[])
{
  bool_t analyze_output, equilibria_only;
  int    niter, nact;

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
  
  readAtomicModels();
  readMolecularModels();
  SortLambda();
  
  Background(analyze_output=TRUE, equilibria_only=FALSE);
  convertScales(&atmos, &geometry);
  getRays(&geometry);
  getBoundary(&geometry);

  getProfiles();
  initSolution();
  initScatter();
  getCPU(1, TIME_POLL, "Total Initialize");

  /* --- Solve radiative transfer for active ingredients -- --------- */

  Iterate(input.NmaxIter, input.iterLimit);

  niter = 0;
  while (niter < input.NmaxScatter) {  
    if (solveSpectrum(FALSE, FALSE) <= input.iterLimit) break;
    niter++;
  }
  /* --- Write output files --                     ------------------ */

  getCPU(1, TIME_START, NULL);
 
  writeInput();
  writeAtmos(&atmos);
  writeGeometry(&geometry);
  writeSpectrum(&spectrum);
  writeFlux("flux.out");

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

  getCPU(1, TIME_POLL, "Output");
  printTotalCPU();
}
/* ------- end ---------------------------- main.c ------------------ */
