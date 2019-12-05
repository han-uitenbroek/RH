/* ---------------------------------------- constant.h --------------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Wed Nov 17 16:27:27 2010 --

       --------------------------                      ----------RH-- */

#ifndef __CONSTANT_H__
#define __CONSTANT_H__

/* --- Defines physical constants.
       The code uses SI definitions throughout the program.

 Note: Wavelengths are stored in nm rather than m.

 Note: When the electron charge is given in Gaussian units with 
       e = 4.803E-10 e.s.u [electrostatic units] the MKSA
       system requires it to be given as e/sqrt(4*PI*EPSILON_0) with e
       in Coulomb.

 Note: \mu_0 / 4\pi \equiv 1.0E-7 

  Ref: http://physics.nist.gov/PhysRefData --          -------------- */

/* --- Physical constants --                           -------------- */

#define  CLIGHT      2.99792458E+08      /* Speed of light [m/s]      */
#define  HPLANCK     6.6260755E-34       /* Planck's constant [Js]    */
#define  KBOLTZMANN  1.380658E-23        /* Boltzman's constant [J/K] */
#define  AMU         1.6605402E-27       /* Atomic mass unit [kg]     */
#define  M_ELECTRON  9.1093897E-31       /* Electron mass [kg]        */
#define  Q_ELECTRON  1.60217733E-19      /* Electron charge [C]       */
#define  EPSILON_0   8.854187817E-12     /* Vacuum permittivity [F/m] */
#define  MU_0        1.2566370614E-06    /* Magnetic induct. of vac.  */
#define  RBOHR       5.29177349E-11      /* Bohr radius [m]           */
#define  E_RYDBERG   2.1798741E-18       /* Ion. pot. Hydrogen [J]    */
#define  EV          1.60217733E-19      /* One electronVolt [J]      */
#define  THETA0      5.03974756E+03      /* log10(e) * eV/k [K^-1]    */
#define  ABARH       7.42E-41            /* polarizabilty of Hydrogen
                                            in [Fm^2]                 */


/* --- Unit conversions --                             -------------- */

#define  NM_TO_M         1.0E-09
#define  CM_TO_M         1.0E-02
#define  KM_TO_M         1.0E+03
#define  ERG_TO_JOULE    1.0E-07
#define  G_TO_KG         1.0E-03
#define  MICRON_TO_NM    1.0E+03
#define  MEGABARN_TO_M2  1.0E-22
#define  ATM_TO_PA       1.0135E+05      /* Atm to Pascal (N/m^2)     */


/* --- Mathematical constants --                       -------------- */

#ifndef PI
#define  PI          3.14159265358979
#endif
#define  SQRTPI      1.77245385090551


/* --- 1/(2sqrt(2)), needed for anisotropy of radiation -- ---------- */

#define  TWOSQRTTWO  0.35355339059327

#define LARMOR (Q_ELECTRON / (4.0*PI*M_ELECTRON)) * NM_TO_M


/* --- Ionization energy Hmin in [J] --                -------------- */

#define E_ION_HMIN  0.754 * EV

#endif /* !__CONSTANT_H__ */

/* ------- end ---------------------------- constant.h -------------- */
