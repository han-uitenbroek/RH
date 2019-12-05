/* ------- file: -------------------------- spectrum.h --------------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Tue Jul  7 14:50:51 2009 --

       --------------------------                      ----------RH-- */

#ifndef __SPECTRUM_H__
#define __SPECTRUM_H__


/* --- Defines structure containing the overall set of wavelengths,
       and the transitions that are active at each of the wavelengths.
       --                                              -------------- */

#define N_MAX_OVERLAP  100

#define VACUUM_TO_AIR_LIMIT  200.0000
#define AIR_TO_VACUUM_LIMIT  199.9352

/* --- File for storage of intensities in case of angle-dependent PRD */

#define IMU_FILENAME  "Imu.dat"


/* --- Stores set of active transitions at each wavelength and pointers
       to temporary storage for opacity and emissivity components --  */

typedef struct{
  int     *Nactiveatomrt, *Nactivemolrt,
          *Nlower, *Nupper, **lower_levels, **upper_levels;
  double  *chi, *eta, *chip, *chi_c, *eta_c, *sca_c, *chip_c;
  AtomicTransition **art;
  MolTransition **mrt;
} ActiveSet;

/* --- Stores emergent intensities and array of active sets -- ------ */

typedef struct {
  bool_t   vacuum_to_air, updateJ;
  int      Nspect, *PRDindex, fd_J, fd_J20, fd_Imu;
  double  *lambda, **J, **I, **Stokes_Q, **Stokes_U, **Stokes_V, **J20;
  ActiveSet *as;
} Spectrum;


/* --- Associated function prototypes --               -------------- */

double Formal(int nspect, bool_t eval_operator, bool_t redistribute);
double solveSpectrum(bool_t eval_operator, bool_t redistribute);

void   addtoGamma(int nspect, double wmu, double *P, double *Psi);
void   addtoRates(int nspect, int mu, bool_t to_obs, double wmu,
		  double *I, bool_t redistribute);
void   initScatter(void);

void   StokesK(int nspect, int k, double chi_I, double K[4][4]);
void   addtoCoupling(int nspect);


/* --- What type of lines are present in active set as? -- ---------- */

bool_t containsPolarized(ActiveSet *as);
bool_t containsBoundBound(ActiveSet *as);
bool_t containsPRDline(ActiveSet *as);
bool_t containsActive(ActiveSet *as);

void init_as(ActiveSet *as);
void alloc_as(int nspect, bool_t crosscoupling);
void free_as(int nspect, bool_t crosscoupling);
void initSpectrum(void);
void Opacity(int nspect, int mu, bool_t top_to_bottom, bool_t activate);
void Planck(int Nspace, double *T, double lambda0, double *Bnu);

void readJlambda(int nspect, double *J);
void writeJlambda(int nspect, double *J);
void readJ20lambda(int nspect, double *J20);
void writeJ20lambda(int nspect, double *J20);

void readImu(int nspect, int mu, bool_t to_obs, double *I);
void writeImu(int nspect, int mu, bool_t to_obs, double *I);
void writeSpectrum(Spectrum *spectrum);
void writeOpacity();


/* --- Wavelength conversion --                         ------------- */

void vacuum_to_air(int Nlambda, double *lambda_vac, double *lambda_air);
void air_to_vacuum(int Nlambda, double *lambda_air, double *lambda_vac);

#endif /* !__SPECTRUM_H__ */

/* ------- end ---------------------------- spectrum.h -------------- */
