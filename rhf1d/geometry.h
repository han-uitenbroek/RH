/* ------- file: -------------------------- geometry.h --------------

       Version:       rh2.0, 1-D plane-parallel
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Thu May 31 14:45:26 2018 --

       --------------------------                      ----------RH-- */

#ifndef __GEOMETRY_H__
#define __GEOMETRY_H__

/* --- Define geometric quantities for 1-D plane-parallel version --  */


enum boundcond  {IRRADIATED, ZERO, THERMALIZED, REFLECTIVE};
enum mass_scale {GEOMETRIC, COLUMN_MASS, TAU500};
enum vertical   {TOP=0, BOTTOM};

typedef struct {
  enum     mass_scale  scale;
  enum     boundcond vboundary[2]; 
  int      Ndep, Nrays;
  double  *height, *cmass, *tau_ref, *mux, *muy, *muz, *wmu, *vel,
         **Itop, **Ibottom;
} Geometry;

/* --- Associated function prototypes --               -------------- */

void convertScales(Atmosphere *atmos, Geometry *geometry);
void getAngleQuad(Geometry *geometry);
void getBoundary(Geometry *geometry);
void MULTIatmos(Atmosphere *atmos, Geometry *geometry);
void writeGeometry(Geometry *geometry);


/* --- Formal solution related --                      -------------- */

double Feautrier(int nspect, int mu, double *chi, double *S,
		 enum FeautrierOrder order, double *P, double *Psi);

void Piecewise_1D(int nspect, int mu, bool_t to_obs,
		  double *chi, double *S, double *I, double *Psi);
void Piecewise_Linear_1D(int nspect, int mu, bool_t to_obs,
			 double *chi, double *S, double *I, double *Psi);

void Piecewise_Bezier_1D(int nspect, int mu, bool_t to_obs,
			 double *chi, double *S, double *I, double *Psi);

/* --- Prototype for parabolic DELO solver --          -------------- */

void Piece_Stokes_1D(int nspect, int mu, bool_t to_obs,
		     double *chi_I, double **S, double **I,
		     double *Psi);


/* --- Prototypes for cubic Bezier solvers --          -------------- */

void Piece_Stokes_Bezier3_1D(int nspect, int mu, bool_t to_obs,
			     double *chi, double **S, double **I,
			     double *Psi);

void Piecewise_Bezier3_1D(int nspect, int mu, bool_t to_obs,
			  double *chi, double *S, double *I,
			  double *Psi);

#endif /* !__GEOMETRY_H__ */

/* ---------------------------------------- geometry.h -------------- */
