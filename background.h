/* ------- file: -------------------------- background.h ------------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Wed Jan 16 14:10:02 2013 --

       --------------------------                      ----------RH-- */

#ifndef __BACKGROUND_H__
#define __BACKGROUND_H__

/* --- Include file for background opacities and emissivities. -- --- */


/* --- Maximum number of iterations and minimum accuracy for chemical
       equilibrium calculations --                     -------------- */

#define N_MAX_CHEM_ITER    10
#define CHEM_ITER_LIMIT    1.0E-3

/* --- Maximum number of iterations and minimum accuracy for hydrostatic
       equilibrium calculations --                     -------------- */

#define N_MAX_HSE_ITER    20
#define HSE_ITER_LIMIT    1.0E-2


/* --- Structure definitions --                        -------------- */

struct Linelist {
  bool_t  used;
  double *adamp;
  AtomicLine *line;
};

/* --- Associated function prototypes --               -------------- */

void   ChemicalEquilibrium(int NmaxIter, double iterLimit);

double Gaunt_bf(double lambda, double n_eff, int charge);
double Gaunt_ff(double lambda, int charge, double T);
void   FMetals(double *Fmetal);
void   dFMetals(double *dFdne);
void   Hydrostatic(int NmaxIter, double iterLimit);

void   Background(bool_t analyze_output, bool_t equilibria_only);
void   backgrOpac(int Nlambda, double *lambda);
bool_t duplicateLevel(Atom *atom, char *label);
bool_t duplicateLine(Atom *atom, char *labeli, char *labelj);
void   readBackground(int nspect, int mu, bool_t top_to_bottom);
int    writeBackground(int nspect, int mu, bool_t to_obs,
		       double *chi_c, double *eta_c, double *sca_c,
		       double *chip_c);
void   writeBRS(void);
void   readBRS(void);

void   readKuruczLines(char *fileName);
int    rlk_ascend(const void *v1, const void *v2);
void   rlk_locate(int N, RLK_Line *lines, double lambda, int *low);

bool_t Hminus_bf(double lambda, double *chi, double *eta);
bool_t Hminus_ff(double lambda, double *chi);
bool_t Hminus_ff_long(double lambda, double *chi);
bool_t Hydrogen_bf(double lambda, double *chi, double *eta);
void   Hydrogen_ff(double lambda, double *chi);
bool_t H2minus_ff(double lambda, double *chi);
bool_t H2plus_ff(double lambda, double *chi);
bool_t Rayleigh(double lambda, Atom *atom, double *chi);
bool_t Rayleigh_H2(double lambda, double *chi);
void   Thomson(double *chi);
bool_t Metal_bf(double lambda, int Nmetal, Atom *metals,
		double *chi, double *eta);
bool_t OH_bf_opac(double lambda, double *chi, double *eta);
bool_t CH_bf_opac(double lambda, double *chi, double *eta);


#endif /* !__BACKGROUND_H__ */

/* ------- end ---------------------------- background.h ------------ */
