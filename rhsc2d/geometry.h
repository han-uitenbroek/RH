/* ------- file: -------------------------- geometry.h --------------

       Version:       rh2.0, 2-D Cartesian
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Dec 17 1996

       --------------------------                      ----------RH-- */

#ifndef __GEOMETRY_H__
#define __GEOMETRY_H__

/* --- Defines structures to contain the angle quadrature and
       interpolation mesh for 2-D Cartesian version. -- ------------- */


enum boundary   {FIXED, PERIODIC};
enum boundval   {IRRADIATED, ZERO, THERMALIZED, REFLECTIVE};
enum horizontal {LEFT=0, RIGHT};
enum vertical   {TOP=0, BOTTOM};
enum direction  {UPWIND=0, DOWNWIND};
enum intersect  {VERTICAL=0, HORIZONTAL};
enum SC_order   {SC_LINEAR, SC_QUADRATIC};
enum triplet    {LOWER_TRIPLET=0, UPPER_TRIPLET};
enum sweep      {DOWN=0, UP};


/* --- The Stencil structure holds the interpolation coefficients and
       indices in the UPWIND and DOWNWIND direction for each ray at
       each spatial location in the grid --            -------------- */ 

typedef struct {
  enum   intersect intersect[2];
  enum   SC_order order[2];
  enum   triplet triplet[2];
  int    index[2][3];
  double coeff[2][3], ds[2], fraction[2];
} Stencil;

typedef struct {
  int      Nlc, lstart;
  Stencil *stencil;
} LongChar;

typedef struct {
  enum     boundary hboundary;
  enum     boundval bvalue[2]; 
  int      Nx, Nz, Nrays;
  double  *x, *dx, *z, *dz, *mux, *muy, *muz, *wmu, *vx, *vz, **Itop,
         **Ibottom, **Ileft, **Iright;
  Stencil  **stencil;
  LongChar **longchar;
} Geometry;


/* --- Associated function prototypes --               -------------- */

void   getAngleQuadr(Geometry *geometry);
void   fillMesh(Geometry *geometry);
bool_t writeFlux(char *fileName);
void   writeGeometry(Geometry *geometry);
void   readAtmos(Atmosphere *atmos, Geometry *geometry);
void   getBoundary(Atmosphere *atmos, Geometry *geometry);
double Quadr(double *values, Stencil *stencil,
	     enum direction direction, bool_t monotonic);
void   StokesK_2D(int nspect, Stencil *st, enum direction direction,
		  bool_t monotonic, double chi_I, double K[4][4]);


/* --- Formal solution related --                      -------------- */


double SolveLong(LongChar *lc, int local, double *chi, double *S,
		 double *I);

double Feautrier(int nspect, int lx, double *chi, double *S,
                 enum FeautrierOrder order, double *P, double *Psi);


void   Piecewise_2D(Geometry *geometry, int nspect, int mu,
		    bool_t to_observer, double *chi, double *S,
		    double *I, double *Psi);
void   Piecewise_Linear_2D(Geometry *geometry, int nspect, int mu,
			   bool_t to_observer, double *chi, double *S,
			   double *I, double *Psi);
void   Piecewise_Bezier_2D(Geometry *geometry, int nspect, int mu,
			   bool_t to_observer, double *chi, double *S,
			   double *I, double *Psi);
void   Piecewise_Hermite_2D(Geometry *geometry, int nspect, int mu,
			    bool_t to_observer, double *chi, double *S,
			    double *I, double *Psi);


/* --- Prototype for parabolic DELO solver --          -------------- */


void   SolveLongStokes(int nspect, LongChar *lc, int local, double *chi,
		       double **S, double **I, double *I_uw);
void   Piece_Stokes_2D(Geometry *geometry,
		       int nspect, int mu, bool_t to_obs, double *chi,
		       double **S_pol, double **I_pol, double *Psi);



/* --- Prototypes for cubic Bezier solvers --          -------------- */

void   Piecewise_Bezier3_2D(Geometry *geometry, int nspect, int mu,
			    bool_t to_observer, double *chi, double *S,
			    double *I, double *Psi);

void   Piece_Stokes_Bezier3_2D(Geometry *geometry, int nspect, int mu,
			       bool_t to_observer, double *chi, double **S,
			       double **I, double *Psi);


#endif /* !__GEOMETRY_H__ */

/* ---------------------------------------- geometry.h -------------- */
