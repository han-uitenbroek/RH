/* ------- file: -------------------------- atom.h ------------------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Mon Sep 22 15:53:02 2025 --

       --------------------------                      ----------RH-- */

#ifndef __ATOM_H__
#define __ATOM_H__


/* --- Defines atomic data structures. --              -------------- */


#define ATOM_LABEL_WIDTH   20
#define ATOM_ID_WIDTH       2
#define MOLECULE_ID_WIDTH  10

#define PRD_QCORE   2.0
#define PRD_QWING   4.0
#define PRD_QSPREAD 5.0
#define PRD_DQ      0.25


enum type        {ATOMIC_LINE, ATOMIC_CONTINUUM,
		  VIBRATION_ROTATION, MOLECULAR_ELECTRONIC};
enum ftype       {FIXED_LINE, FIXED_CONTINUUM};
enum Trad_option {TRAD_ATMOSPHERIC, TRAD_PHOTOSPHERIC, TRAD_CHROMOSPHERIC};
enum vdWaals     {UNSOLD, RIDDER_RENSBERGEN, BARKLEM, KURUCZ};
enum fit_type    {KURUCZ_70, KURUCZ_85, SAUVAL_TATUM_84, IRWIN_81, TSUJI_73};
enum Hund        {CASE_A, CASE_B};
enum Barklemtype {SP, PD, DF};
enum orbit_am    {S_ORBIT=0, P_ORBIT, D_ORBIT, F_ORBIT};
enum zeeman_cpl  {LS_COUPLING=0, JK_COUPLING, JJ_COUPLING};

/* --- Structure prototypes --                         -------------- */

typedef struct Atom Atom;
typedef struct Molecule Molecule;
typedef struct AtomicTransition AtomicTransition;
typedef struct MolTransition MolTransition;
typedef struct AtomicLine AtomicLine;
typedef struct ZeemanMultiplet ZeemanMultiplet;
typedef struct rhthread rhthread;
typedef struct Paschenstruct Paschenstruct;

/* --- Structure defines radiative transition --       -------------- */

struct AtomicLine {
  bool_t   symmetric, polarizable, Voigt, PRD;
  enum vdWaals vdWaals;
  int      i, j, Nlambda, Nblue, Ncomponent, Nxrd, fd_profile,
           **id0, **id1;
  double   lambda0, *lambda, isotope_frac, g_Lande_eff,
           Aji, Bji, Bij, *Rij, *Rji, **phi, **phi_Q, **phi_U, **phi_V,
         **psi_Q, **psi_U, **psi_V, *wphi, *Qelast, Grad, cvdWaals[4],
    cStark, qcore, qwing, **rho_prd, *c_shift, *c_fraction, **frac;
  FILE    *fp_GII;
  struct Ng *Ng_prd;
  Atom *atom;
  AtomicLine **xrd;
  pthread_mutex_t rate_lock;
  ZeemanMultiplet *zm;
};

typedef struct {
  bool_t  hydrogenic;
  int     i, j, Nlambda, Nblue;
  double  lambda0, *lambda, isotope_frac, alpha0, *alpha, *Rij, *Rji;
  Atom *atom;
  pthread_mutex_t rate_lock;
} AtomicContinuum;

typedef struct {
  enum type type;
  enum Hund Hundi, Hundj;
  bool_t symmetric, Voigt, polarizable;
  char   configi[3], configj[3], parityi[2], parityj[2];
  int    vi, vj, Nlambda, Nblue, subi, subj, Lambdai, Lambdaj, ecnoi, ecnoj;
  double lambda0, *lambda, isotope_frac, Ei, Ej, gi, gj, Si, Sj,
         Omegai, Omegaj, **phi, *wphi, g_Lande_eff,
         Grad, qcore, qwing, Aji, Bji, Bij;
  Molecule *molecule;
  ZeemanMultiplet *zm;
} MolecularLine;

typedef struct {
  enum   ftype type;
  enum   Trad_option option;
  int    i, j;
  double lambda0, strength, Trad;
  Atom *atom;
} FixedTransition;

struct AtomicTransition {
  enum type type;
  union {
    AtomicLine *line;
    AtomicContinuum *continuum;
  } ptype;
  Atom *atom;
};

struct MolTransition {
  enum type type;
  union {
    MolecularLine *vrline;
  } ptype;
  Molecule *molecule;
};

struct rhthread {
  double **gij, **Vij, **wla, **chi_up, **chi_down, **Uji_down, *eta;
};

struct Atom {
  char    ID[ATOM_ID_WIDTH+1], **label, *popsinFile, *popsoutFile;
  bool_t  active, NLTEpops;
  enum solution initial_solution;
  int     Nlevel, Nline, Ncont, Nfixed, Nprd, *stage, periodic_table,
          activeindex;
  long    offset_coll;
  double  abundance, weight, *g, *E, **C, *vbroad, **n, **nstar,
         *ntotal, **Gamma;  
  AtomicLine *line;
  AtomicContinuum *continuum;
  FixedTransition *ft;
  struct Ng *Ng_n;
  rhthread *rhth;
  pthread_mutex_t Gamma_lock;
  FILE *fp_input;
};

typedef struct {
  char     ID[ATOM_ID_WIDTH+1];
  bool_t   abundance_set;
  int     *mol_index, Nstage, Nmolecule;
  double   weight, abund, *ionpot, **pf, **n;
  Atom *model;
} Element;

struct Molecule {
  char    ID[MOLECULE_ID_WIDTH+1], *popsFile, *configs;
  bool_t  active;
  enum    fit_type fit;
  enum    solution initial_solution;
  int    *pt_index, *pt_count, Nelement, Nnuclei, Npf, Neqc, Nrt,
          charge, Nconfig, Nv, NJ, activeindex;
  double  Ediss, Tmin, Tmax, weight, *vbroad, *pf_coef, *eqc_coef,
         *pf, **pfv, *n, **nv, **nvstar, *C_ul, **Gamma;
  MolecularLine *mrt;
  struct Ng *Ng_nv;
  rhthread *rhth;
  pthread_mutex_t Gamma_lock;
};

typedef struct {
  int    L, L1, l1, l2, l;
  double g, E, S, J, S1, J1, j1, j2, K, gL, hfs;
  enum zeeman_cpl cpl;
  bool_t zm_explicit;
} RLK_level;
  
typedef struct {
  bool_t polarizable;
  enum vdWaals vdwaals;
  int    pt_index, stage, isotope;
  double lambda0, Bji, Aji, Bij,
         Grad, GStark, GvdWaals, hyperfine_frac,
         isotope_frac, iso_dl, cross, alpha;
  RLK_level level_i, level_j;
  ZeemanMultiplet *zm;
} RLK_Line;

struct ZeemanMultiplet{
  int     Ncomponent, *q;
  double *shift, *strength, g_eff;
};

typedef struct {
  int     N1, N2;
  double *neff1, *neff2, **cross, **alpha;
} Barklemstruct;

struct Paschenstruct{
  int Nj;
  double *eigenval;
  double **C;
};


/* --- Associated function prototypes --               -------------- */

void   initSolution(void);
void   Iterate(int NmaxIter, double iterLimit);

void   readAtomicModels(void);
void   readMolecularModels(void);
void   statEquil(Atom *atom, int isum);
double updatePopulations(int niter);

void CollisionRate(Atom *atom, FILE *atomFile);
void Damping(AtomicLine *line, double *adamp);
void FixedRate(Atom *atom);
void freeAtom(Atom *atom);
void freeAtomicLine(AtomicLine *line);
void freeAtomicContinuum(AtomicContinuum *continuum);
void getfjk(Element *element, double ne, int k, double *fjk, double *dfjk);
void getLambda(AtomicLine *line);
void getLambdaCont(AtomicContinuum *continuum, double lambdamin);
double getwlambda_line(AtomicLine *line, int la);
double getwlambda_cont(AtomicContinuum *continuum, int la);
void initAtom(Atom *atom);
void initAtomicLine(AtomicLine *line);
void initAtomicContinuum(AtomicContinuum *continuum);

void initZeeman(ZeemanMultiplet *zm);
void freeZeeman(ZeemanMultiplet *zm);
double zm_gamma(double J, double S, double L);

void initGammaAtom(Atom *atom, int iter);
void initGammaMolecule(Molecule *molecule);

void LTEpops(Atom *atom, bool_t Debeye);
void LTEpops_elem(Element *element);
void SetLTEQuantities(void);
void getProfiles(void);
void Profile(AtomicLine *line);
void readProfile(AtomicLine *line, int lamu, double *phi);
void writeProfile(AtomicLine *line, int lamu, double *phi);
void readAtom(Atom *atom, char *atomFileName);
void readPopulations(Atom *atom);
void SortLambda(void);
void Stark(AtomicLine *line, double *GStark);
void StarkLinear(AtomicLine *line, double *GStark);
void VanderWaals(AtomicLine *line, double *GvdW);
void writeAtom(Atom *atom);
void writePopulations(Atom *atom);
void zeroRates(bool_t redistribute);

bool_t readRadRate(Atom *atom);
bool_t writeRadRate(Atom *atom);
bool_t readCollisionRate(Atom *atom);
bool_t writeCollisionRate(Atom *atom);
bool_t readDamping(Atom *atom);
bool_t writeDamping(Atom *atom);

bool_t readBarklemTable(enum Barklemtype type, Barklemstruct *bs);
bool_t getBarklemcross(Barklemstruct *bs, RLK_Line *rlk);
bool_t getBarklemactivecross(AtomicLine *line);
bool_t getBarklemExplicit(AtomicLine *line);
bool_t determinate_abo(char *label,  int *l);


/* --- Associated function prototypes --               -------------- */

void statEquilMolecule(Molecule *molecule, int isum);

void   COcollisions(Molecule *molecule);
void   H2collisions(Molecule *molecule);
void   MolecularDamping(MolecularLine *mrt, double *adamp);
double equilconstant(Molecule *molecule, double T);
void   freeMolecule(Molecule *molecule);
void   freeMolLine(MolecularLine *mrt);
double getwlambda_mrt(MolecularLine *mrt, int la);
void   initMolecule(Molecule *molecule);
void   initMolLine(MolecularLine *mrt, enum type line_type);
void   MolecularProfile(MolecularLine *mrt);
void   mrt_locate(int N, MolecularLine *lines, double lambda, int *low);
void   LTEmolecule(Molecule *molecule);
double partfunction(Molecule *molecule, double T);
void   readMolecule(Molecule *molecule, char *fileName, bool_t active);
void   readMolecularLines(Molecule *molecule, char *line_data);
void   readMolPops(Molecule *molecule);
void   writeMolPops(Molecule *molecule);


/* --- Redistribution function --                      -------------- */

double GII(double adamp, double waveRatio, double q_emit, double q_abs);
double PII(double adamp, double waveRatio, double q_emit, double q_abs);
double RII(double v_emit, double v_abs, double adamp, int mu1, int mu2);
void   Redistribute(int NmaxIter, double iterLimit);
void   PRDScatter(AtomicLine *line, enum Interpolation representation);
void   PRDAngleScatter(AtomicLine *PRDline,
		       enum Interpolation representation);
void   PRDAngleApproxScatter(AtomicLine *PRDline,
			     enum Interpolation representation);


/* --- Polarization related --                         -------------- */

void   adjustStokesMode();
bool_t determinate(char *label, double g, int *n, double *S, int *L,
		   double *J);
double effectiveLande(AtomicLine *line);
double Lande(double S, int L, double J);
void   StokesProfile(AtomicLine *line);
void   Zeeman(AtomicLine *line);
void   MolZeeman(MolecularLine *mrt);
double MolLande_eff(MolecularLine *mrt);
int    getOrbital(char orbit);
double ZeemanStrength(double Ju, double Mu, double Jl, double Ml);

double w3js(double j1, double j2, double j3,
	    double m1, double m2, double m3);
double w6js(double j1, double j2, double j3,
	    double J1, double J2, double J3);


#endif /* !__ATOM_H__ */

/* ------- end ---------------------------- atom.h ------------------ */
