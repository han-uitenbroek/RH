/* ------- file: -------------------------- geometry.h --------------

       Version:       rh2.0, 1-D spherically symmetric
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Fri Jan 28 11:02:03 2005 --

       --------------------------                      ----------RH-- */

#ifndef __GEOMETRY_H__
#define __GEOMETRY_H__

/* --- Define Geometry structure and prototypes for spherical
       geometry --                                     -------------- */


#define  PTOP  2

enum boundcond     {ZERO, THERMALIZED, IRRADIATED};
enum depthscale    {GEOMETRIC, COLUMN_MASS, TAU500};
enum location      {TOP=0, CORE};
enum raytype       {SHELL_RAY, CORE_RAY};

/* --- The Ray structure has the following components:

         type: enum of type raytype indicating whether this ray is
               a grazing (type = SHELL_RAY) or core ray (type = CORE_RAY)
           Ns: Number of grid points along ray
	    p: Impact parameter of ray (absolute value from center, in m) 
            s: Pathlength along ray (counted from the centerline)
          xmu: Array of direction cosines along ray.
          wmu: Array of integration weights along ray.
	  --                                           -------------- */

typedef struct {
  enum   raytype type;
  int    Ns;
  double p, *s, *xmu, *wmu;
} Ray;

typedef struct {
  enum    depthscale scale;
  enum    boundcond rboundary[2];
  int     Nradius, Nrays, Ncore, Ninter;
  double  Radius, *r, *cmass, *tau_ref, *vel, *Itop;
  Ray *rays;
} Geometry;

/* --- Associated function prototypes --               -------------- */

void   addtoGamma_sphere(int nspect, Ray *ray, double *P, double *Psi);
void   addtoRates_sphere(int nspect, Ray *ray, double *P,
			 bool_t redistribute);
void   convertScales(Atmosphere *atmos, Geometry *geometry);
double Feautrier(int nspect, int mu, double *chi, double *S,
		 enum FeautrierOrder F_order, double *P, double *Psi);
void   getRays(Geometry *geometry);
void   getBoundary(Geometry *geometry);
void   readAtmos(Atmosphere *atmos, Geometry *geometry);
bool_t writeFlux(char *fileName);
void   writeGeometry(Geometry *geometry);

#endif /* !__GEOMETRY_H__ */

/* ---------------------------------------- geometry.h -------------- */
