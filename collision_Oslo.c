/* ------- file: -------------------------- collision.c -------------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Sat Sep 19 15:54:07 2009 --

       --------------------------                      ----------RH-- */

/* --- Reads collisional data from atomic data file and computes
       collisional rates.

       Format is similar to MULTI's routine GENCOL:

         Input line can be either

           TEMP  Nitem  T[0]     ...   T[Nitem-1]

         to read in the temperature grid for collisional coefficients, or

           KEYWORD  i1  i2   coeff[0]    ...   coeff[Nitem-1]

         to read the coefficients for transitions between i1 and i2.
         Multiple entries of the temperature grid are allowed, but
         at least one entry with the correct number of grid points
         has to precede an entry with coefficients.

         Allowed keywords are:

         Keyword    Transition type
         ----------------------------------------------------

         TEMP  -->  Temperature grid

         OMEGA -->  Collisional de-excitation of ions by electrons
         CE    -->  Collisional de-excitation of neutrals by electrons
         CI    -->  Collisional ionization by electrons
         CP    -->  Collisional de-excitation by protons

         CH0   -->  Charge exchange of ion with neutral hydrogen
         CH+   -->  Charge exchange of neutral with protons

         END   -->  End of input data
         ----------------------------------------------------

 Note: Unit of number density is m^-3.

       Convention: C_ij = C[i][j] represents the
                   transition j --> i

       --                                              -------------- */

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "constant.h"
#include "error.h"
#include "inputs.h"
#include "statistics.h"

#define COMMENT_CHAR "#"
#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )

/* --- Function prototypes --                          -------------- */
double E1(double x);
double fone(double x);
double ftwo(double x);
void atomnm(int anr,char *cseq);
int atomnr(char id[ATOM_ID_WIDTH+1]);
void ar85cea(int i,int j,int k, double *cea,struct Atom *atom);
void rowcol(int i, int *row, int *col);

/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern char messageStr[];


/* ------- begin --------------------------------- rowcol.c --------- */
void rowcol(int anr, int *row, int *col){
/*  */
  int istart[7]={0,2,10,18,36,54,86};
  int i;

  for (i=0; i<5; i++) {
    if (anr > istart[i] && anr <= istart[i+1]) {
      *row=i;
      break;
    }
  }
  *col=anr-istart[*row];
  *row = *row + 1;
  
 }
/* ------- end ----------------------------------- rowcol.c --------- */



/* ------- begin ----------------------------------- fone.c --------- */
double fone(double x){
/* function f_1 from Arnoud & Rothenflug, 1985, A&ASS, 60,425  */

  double y;

  /*  if (x<=0.02) {
    y=exp(x) * (-1*log(x) -0.5772 + x);
  } else if (x<=1.5) {
    y = log((x+1.0)/x) + (0.36+0.03*pow(x+0.01,-0.5)) / pow(x-1,2);
  } else if  (x<=10.0) {
    y = log((x+1.0)/x) + (0.36+0.03*pow(x+0.01,0.5)) / pow(x-1,2);
  } else if  (x<=50.0) {
    y= 1.0/x * (1.0+1.0/x +2.0/pow(x,2) +6.0/pow(x,3) +24.0/pow(x,4));
  } else {
      y = 1.0/x;
      } */
  if (x<=50.0) { 
    y = exp(x)*E1(x);
  } else{
    y= 1.0/x;
  }


 return y;

}
/* ------- end ------------------------------------- fone.c --------- */

/* ------- begin ----------------------------------- ftwo.c --------- */
double ftwo(double x){
/* function f_2 from Arnoud & Rothenflug, 1985, A&ASS, 60,425, with
   improved description when x<4 from Hummer, 1983, jqsrt, 30 281  */

  double y;
  double brk=4.0;

  double p[15]={1.0000e+00,2.1658e+02,2.0336e+04,1.0911e+06,3.7114e+07,
		8.3963e+08,1.2889e+10,1.3449e+11,9.4002e+11,4.2571e+12,
		1.1743e+13,1.7549e+13,1.0806e+13,4.9776e+11,0.0000};
  double q[15]={1.0000e+00,2.1958e+02,2.0984e+04,1.1517e+06,4.0349e+07,
		9.4900e+08,1.5345e+10,1.7182e+11,1.3249e+12,6.9071e+12,
		2.3531e+13,4.9432e+13,5.7760e+13,3.0225e+13, 3.3641e+12};

  int i;

  double px,xfact,qx,gamma,f0x,count,fact,term;

  if (x>brk) {

    px=p[0];
    xfact=1.0;
    for (i=1;  i<15;  i++) {
      xfact = xfact/x;
      px=px+p[i]*xfact;
    }
     
    qx=q[0];
    xfact=1.0;
    for (i = 1;  i<15;  i++) {
      xfact = xfact/x;
      qx=qx+ q[i]*xfact;
    }

    y=px/qx/x/x;

  } else {

    gamma=0.5772156649;
    f0x = PI*PI/12.;
    term=1.0;
    count=0.0;
    fact=1.0;
    xfact=1.0;
    while (fabs(term/f0x) > 1.e-8){
      count=count+1.0;
      fact = fact*count;
      xfact=xfact*(-x);
      term = xfact/count/count/fact;
      f0x=f0x+term;
      if (count>100.) /*ERROR CODE HERE */ ;
    }
    y=exp(x)*((log(x)+gamma)*(log(x)+gamma)*0.5 + f0x);
  }

 return y;

}
/* ------- end ------------------------------------- ftwo.c --------- */



/* ------- begin --------------------------------- atomnr.c --------- */
void atomnm(int anr,char *cseq){
  /*  atomnm 94-02-22  new routine: (philip judge)
      gives atomic name of arelement if i is an integer containing
      e.g. if input  is i it will return 'h ', 2, it will return 'he', etc. */  
  
  int i;
  
  char *arelem[29] =
    {"H ","HE","LI","BE","B ","C ","N ","O ","F ","NE", 
     "NA","MG","AL","SI","P ","S ","CL","AR", 
     "K ","CA","SC","TI","V ","CR","MN","FE","CO","NI","  "};
  
  i=anr-1;
  
  if (i<=0 || i<28) {
    strcpy(cseq,arelem[i]);
  } else {
    strcpy(cseq,arelem[28]);
  }
   
}
/* ------- end ----------------------------------- atomnm.c --------- */


/* ------- begin -------------------------------- ar85cea.c --------- */
void ar85cea(int i,int j,int k, double *cea,struct Atom *atom) {
  
/*   new routine for computing collisional autoionization 
     rates using formalism and formulae from arnaud and rothenflug 1985
     ar85cea 94-02-22  new routine: (philip judge)
     
     96-03-07  modifications: (philip judge)
     bug fixed: cup initialized to zero    */

  char cseq[ATOM_ID_WIDTH+1];
  int iz,ichrge,isoseq;
  double zz,cup,bkt,b,zeff,iea,y,f1y,a,g;
  
  /*  initialize output to zero  */
  *cea=0.0;
  y=0.0;
  f1y=0.0;
  cup=0.0;

  /*  find element*/
  iz=atomnr(atom->ID);
  zz=iz;
     
  if (iz<1 || iz>92) {
    /* ERROR CODE HERE */
  }

  /*  find iso-electronic sequence */
  isoseq=iz-atom->stage[i];
  atomnm(isoseq,&cseq[0]);
  ichrge=atom->stage[i];
  
  /* temperature in eV */
  bkt=KBOLTZMANN*atmos.T[k]/EV;

 
  /*printf("%i %i %i %i %i %e %s\n",i,j,k,iz,isoseq,zz,cseq);*/

  /* lithium sequence */
  if ( !strcmp(cseq,"LI") ) { 
    
    iea=13.6*(pow(zz-0.835,2) - 0.25*pow(zz-1.62,2));
    b=1. / (1. + 2.e-4*pow(zz,3));
    zeff=zz-0.43;
    y=iea/bkt;
    f1y=fone(y);
    g= 2.22*f1y + 0.67*(1.-y*f1y) + 0.49*y*f1y +1.2*y*(1.-y*f1y);
    
    cup=1.60e-07*1.2*b/pow(zeff,2)/sqrt(bkt)*exp(-y)*g * pow(CM_TO_M,3);
    
    /* special cases */
 
    /* HU -- I think there are a couple of cut & paste errors: C instead
       of N and O. For C this will lead to an erroneous factor of 1.25
       instead of 0.6 */
   
    /* c iv - app a ar85 */
    if (!strcmp(atom->ID,"C")) cup = cup*0.6;
    /* n v  - app a ar85 */
    if (!strcmp(atom->ID,"C")) cup = cup*0.8;
    /* o vi - app a ar85 */
    if (!strcmp(atom->ID,"C")) cup = cup*1.25;  
    
    /* sodium sequence */
  } else if (!strcmp(cseq,"NA") ) {
    
    if (iz <= 16) {
      
      iea=26.*(zz-10.);
      a=2.8e-17*pow(zz-11.,-0.7);
      y=iea/bkt;
      
      f1y=fone(y);
      cup= 6.69e+7 * a * iea / sqrt(bkt) * exp(-y)* (1. - y*f1y) * pow(CM_TO_M,3);
      
     

    } else if  ( iz>=18 && iz<=28 ) {
      
      iea=11.*(zz-10.)*sqrt(zz-10.);	 
      a=1.3e-14*pow(zz-10.,-3.73);
      y=iea/bkt;
      f1y=fone(y);
      cup= 6.69e+7 * a * iea / sqrt(bkt) * exp(-y)  
	*(1. - 0.5*(y - SQ(y) + SQ(y)*y*f1y));
      
    } else { 

      cup=0.0;
      
    }
    
  }
  
  /* magnesium-sulfur sequences */

  if (!strcmp(cseq,"MG") || !strcmp(cseq,"AL") || !strcmp(cseq,"SI") ||
      !strcmp(cseq,"P") || !strcmp(cseq,"S") ) {
    
    if (!strcmp(cseq,"MG")) iea=10.3*pow(zz-10.,1.52);
    if (!strcmp(cseq,"AL")) iea=18.0*pow(zz-11.,1.33);
    if (!strcmp(cseq,"SI")) iea=18.4*pow(zz-12.,1.36);
    if (!strcmp(cseq,"P" )) iea=23.7*pow(zz-13.,1.29);
    if (!strcmp(cseq,"S" )) iea=40.1*pow(zz-14.,1.1 );
    
    a=4.0e-13/SQ(zz)/iea;
    y=iea/bkt;
    f1y=fone(y);
    cup = 6.69e+7 * a * iea / sqrt(bkt) * exp(-y)   
      *( 1.0 - 0.5*(y -SQ(y) + SQ(y)*y*f1y) ) * pow(CM_TO_M,3) ;
  }

  /*  special cases */
  
  if( !strcmp(atom->ID,"CA") && ichrge==0 ) {
    iea= 25.;
    a = 9.8e-17;
    b=1.12;
    cup=6.69e+7 * a * iea / sqrt(bkt) * exp(-y)*(1. + b*f1y)*pow(CM_TO_M,3);
  } else if ( !strcmp(atom->ID,"CA") && ichrge==1 ) {
    a = 6.0e-17;
    iea= 25.;
    b=1.12;
    cup=6.69e+7 * a * iea / sqrt(bkt) * exp(-y)*(1. + b*f1y)*pow(CM_TO_M,3);
  } else if ( !strcmp(atom->ID,"FE") && ichrge==3) {
    a = 1.8e-17;
    iea= 60.;
    b=1.0;
    cup=6.69e+7 * a * iea / sqrt(bkt) * exp(-y)*(1. + b*f1y)*pow(CM_TO_M,3);
  } else if ( !strcmp(atom->ID,"FE") && ichrge==4 ) {
    a = 5.0e-17;
    iea= 73.;
    b=1.0;
    cup=6.69e+7 * a * iea / sqrt(bkt) * exp(-y)*(1. + b*f1y)*pow(CM_TO_M,3);
  }

  *cea=cup;  

}
/* --------- end -------------------------------- ar85cea.c --------- */


/* ------- begin ---------------------------------summers.c --------- */
double summers(int i, int j, double nne, struct Atom *atom){
  /*   
  
  22-Jun-1994 changes begin P.G.Judge ***************************************

  DENSITY SENSITIVE DIELECTRONIC RECOMBINATION

  the term adi may be multiplied by a density-sensitive factor
  if needed- this is crucial for Li and B-like ions colliding with
  impacting electrons.
  This simple formulation was derived from a study of the dependence of
  the dielectronic "bump" in the figures of Summers 1974 
  (Appleton Laboratory internal memo), and fitting according to the
  parameter Ne / z^7
  This should be accurate to typically +/- 0.1 in log in regions
  where it matters.  Worse case is e.g. C like Neon where it underestimates
  density factor by maybe 0.25 in log.

  June 24, 2006 changes begin P.G.Judge 
  original (pre MAY 2006) code
           define rho = nne/ z^7 where z is charge on recombining ion
         rho=10.^(alog10(nne) - 7.* alog10(charge))
         rho0=2.e3
         if(isos eq 'LI' or isos eq 'NA' or isos eq 'K') then rho0 = 3.e1      
         ne_factor = 1./(1. + rho/rho0)^0.14
         print,'ne_factor',ne_factor

 June 24, 2006, more accurate version. 
 get the row and column of the recombined ion's isoelectronic 
 sequence  in the periodic table
 the following parameters mimic the tables 1-19 of
 H. Summers' Appleton Lab Report 367, 1974

Mar 9 2012 Jorrit Leenaarts: stole routine from Phil Judge's DIPER IDL
package

  */


  double y,zz,rho0,rhoq,x,beta;
  int iz,isoseq,row,col;
  char cseq[ATOM_ID_WIDTH+1];
  
  /* find atomic number of element */
  iz= atomnr(atom->ID);
  if (iz<1 || iz>92) {/*ERROR CODE HERE*/}

  /* Charge of recombining ion */
  zz=atom->stage[j];

  /* find iso-electronic sequence of recombined ion */
  isoseq=iz-atom->stage[i];
  atomnm(isoseq,&cseq[0]);

  /* row and column in periodic table */
  rowcol(isoseq,&row,&col);

  rhoq= nne*pow(CM_TO_M,3)/pow(zz,7);
  x= (0.5*zz + (col-1.0)) * row / 3.0;
  beta = -0.2/log(x+2.71828);
  rho0 = 30.0+50.0*x;
  y=pow(1. + rhoq/rho0,beta);

  /*   printf("summers %e %e %e %e \n",nne*pow(CM_TO_M,3), rho0,rhoq,y); */

  return y;

}
/* ------- end ---------------------------------- summers.c --------- */


/* ------- begin --------------------------------- atomnr.c --------- */
int atomnr(char id[ATOM_ID_WIDTH+1]){
  /* Returns atomic number of element with name id. */  
  
  int i,y;
  int n=28;
  char arelem[28][ATOM_ID_WIDTH+1] =
    {"H ","HE","LI","BE","B ","C ","N ","O ","F ","NE", 
     "NA","MG","AL","SI","P ","S ","CL","AR", 
     "K ","CA","SC","TI","V ","CR","MN","FE","CO","NI"};
  
  y=0;
  for (i=0; i<n; i++) {
    if (!strcmp(id,arelem[i])) {
      y=i+1;
      break;
    }
  }
  return y;
}
/* ------- end ----------------------------------- atomnr.c --------- */

/* ------- begin -------------------------- CollisionRate.c --------- */

void CollisionRate(struct Atom *atom, FILE *fp_atom)
{
  const char routineName[] = "CollisionRate";
  register int k, n,m;

  char    inputLine[MAX_LINE_SIZE], *keyword, *pointer,
          labelStr[MAX_LINE_SIZE];
  bool_t  hunt, exit_on_EOF;
  int     nitem, i1, i2, i, j, ij, ji, Nlevel = atom->Nlevel, Nitem,ii,
    status;
  long    Nspace = atmos.Nspace;
  double  dE, C0, *T, *coeff, *C, Cdown, Cup, gij, *np,xj,fac,fxj;
  /* JL additions start */
  double **cdi;
  double **badi;
  int    Ncoef,mshell=5,sumscl;
  double  acolsh,tcolsh,aradsh,xradsh,adish,bdish,t0sh,t1sh,summrs,tg,cdn,cup;
  double  ar85t1,ar85t2,ar85a,ar85b,ar85c,ar85d,t4;
  double de,zz,betab,cbar,dekt,dekti,wlog,wb;
  /* JL additions end */
  getCPU(3, TIME_START, NULL);

  C0 = ((E_RYDBERG/sqrt(M_ELECTRON)) * PI*SQ(RBOHR)) *
    sqrt(8.0/(PI*KBOLTZMANN));

  atom->C = matrix_double(SQ(Nlevel), Nspace);
  for (ij = 0;  ij < SQ(Nlevel);  ij++) {
    for (k = 0;  k < Nspace;  k++) {
      atom->C[ij][k] = 0.0;
    }
  }
  C = (double *) malloc(Nspace * sizeof(double));

  T = coeff = NULL;
  while ((status = getLine(fp_atom, COMMENT_CHAR,
		  inputLine, exit_on_EOF=FALSE)) != EOF) {
    keyword = strtok(inputLine, " ");

    if (!strcmp(keyword, "TEMP")) {

      /* --- Read temperature grid --                  -------------- */

      Nitem = atoi(strtok(NULL, " "));
      T = (double *) realloc(T, Nitem*sizeof(double));
      for (n = 0, nitem = 0;  n < Nitem;  n++) {
        if ((pointer = strtok(NULL, " ")) == NULL) break;
	nitem += sscanf(pointer, "%lf", T+n);
      }
    } else if (!strcmp(keyword, "OMEGA") || !strcmp(keyword, "CE") ||
	       !strcmp(keyword, "CI")    || !strcmp(keyword, "CP") ||
	       !strcmp(keyword, "CH0")   || !strcmp(keyword, "CH+")||
	       !strcmp(keyword, "CH") ) {

      /* --- Read level indices and collision coefficients -- ------- */

      i1 = atoi(strtok(NULL, " "));
      i2 = atoi(strtok(NULL, " "));
      coeff = (double *) realloc(coeff, Nitem*sizeof(double));
      for (n = 0, nitem = 0;  n < Nitem;  n++) {
        if ((pointer = strtok(NULL, " ")) == NULL) break;
	nitem += sscanf(pointer, "%lf", coeff+n);
      }
      /* --- Transitions i -> j are stored at index ji, transitions
	     j -> i are stored under ij. --            -------------- */

      i  = MIN(i1, i2);
      j  = MAX(i1, i2);
      ij = i*Nlevel + j;
      ji = j*Nlevel + i;

    } else if (!strcmp(keyword, "AR85-CHP") || !strcmp(keyword, "AR85-CHH")) {
      
      i1 = atoi(strtok(NULL, " "));
      i2 = atoi(strtok(NULL, " "));
      
      Nitem=6;
      coeff = (double *) realloc(coeff, Nitem*sizeof(double));
      
      for (n = 0, nitem = 0;  n < Nitem;  n++) {
        if ((pointer = strtok(NULL, " ")) == NULL) break;
	nitem += sscanf(pointer, "%lf", coeff+n);
      }

      i  = MIN(i1, i2);
      j  = MAX(i1, i2);
      ij = i*Nlevel + j;
      ji = j*Nlevel + i;

   } else if (!strcmp(keyword,"AR85-CEA") || !strcmp(keyword, "BURGESS")) {

      i1 = atoi(strtok(NULL, " "));
      i2 = atoi(strtok(NULL, " "));      
      
      nitem=1;
      Nitem=1;
      coeff = (double *) realloc(coeff, Nitem*sizeof(double));
      coeff[0] =   atof(strtok(NULL, " "));    
      
      i  = MIN(i1, i2);
      j  = MAX(i1, i2);
      ij = i*Nlevel + j;
      ji = j*Nlevel + i;

    } else if (!strcmp(keyword, "SHULL82")) {
      
      i1 = atoi(strtok(NULL, " "));
      i2 = atoi(strtok(NULL, " "));
      
      Nitem=8;
      coeff = (double *) realloc(coeff, Nitem*sizeof(double));
      
      for (n = 0, nitem = 0;  n < Nitem;  n++) {
        if ((pointer = strtok(NULL, " ")) == NULL) break;
	nitem += sscanf(pointer, "%lf", coeff+n);
      }
      
      i  = MIN(i1, i2);
      j  = MAX(i1, i2);
      ij = i*Nlevel + j;
      ji = j*Nlevel + i;
      
      keyword="SHULL82"; 
      Nitem=8;
      /*BADNELL recipe for dielectronic recombination: bhavna Rathore: 20 Jan 2014 */
    } else if (!strcmp(keyword,"BADNELL")) {
      i1 = atoi(strtok(NULL, " "));
      i2 = atoi(strtok(NULL, " "));
      Ncoef = atoi(strtok(NULL, " "));
      Nitem=Ncoef;
      nitem=Ncoef;
      badi=matrix_double(Ncoef, Nitem);
      for (m = 0; m <= 1;  m++) {
	
        /* HU -- here keyword gets overwritten. This is more a problem of
	   the original code that is exposed here */

	status=getLine(fp_atom, COMMENT_CHAR,inputLine, exit_on_EOF=FALSE);
	badi[m][0]  = atof(strtok(inputLine, " "));
	for (n = 1;  n < Nitem;  n++) {
	  badi[m][n]  = atof(strtok(NULL, " "));
	
	}
      }
        i  = MIN(i1, i2);
      j  = MAX(i1, i2);
      ij = i*Nlevel + j;
      ji = j*Nlevel + i;

      /* HU -- This is very dangerous in C. It assigns a random pointer 
               (pointing to the constant string "BADNELL") to the keyword
               pointer. Memory will get overwritten depending on compiler
               and/or OS */

       keyword="BADNELL"; 

       /* HU -- Nitem is supposed to be the number of items to be read,
                nitem, the number of actually read items. This doesn't
                break it, but doesn't help in catching input errors either */
       Nitem=Ncoef;
       nitem=Ncoef;
      /**************Badnell END */
    } 
      /*switch for density dependent DR coefficent:*/
     else if (!strcmp(keyword, "SUMMERS")) {
      /*  give default multiplication factor of summers density dependence of dielectronic recombination
	  sumscl=0.0 means there is no density dependence
	  sumscl=1.0 means full summers density dependence */
	sumscl = atoi(strtok(NULL, " "));
      }else if (!strcmp(keyword, "AR85-CDI")) {
	
      i1 = atoi(strtok(NULL, " "));
      i2 = atoi(strtok(NULL, " "));
      Ncoef = atoi(strtok(NULL, " "));
      
      if (Ncoef > mshell) {
	sprintf(messageStr, "Ncoef: %i greater than mshell %i",Ncoef, mshell );
	Error(ERROR_LEVEL_2, routineName, messageStr);
      }

      Nitem=5;
      cdi=matrix_double(Ncoef, Nitem);
      
      for (m = 0; m < Ncoef;  m++) {
	
        /* HU -- Same problem with keyword as above */

	status=getLine(fp_atom, COMMENT_CHAR,inputLine, exit_on_EOF=FALSE);
	cdi[m][0]  = atof(strtok(inputLine, " "));
	for (n = 1;  n < Nitem;  n++) {
	  cdi[m][n]  = atof(strtok(NULL, " "));
	}
      }

      i  = MIN(i1, i2);
      j  = MAX(i1, i2);
      ij = i*Nlevel + j;
      ji = j*Nlevel + i;

      /* HU -- Same problems as above with Nitem, nitem and keyword */
      nitem=5;
      Nitem=5;
      keyword="AR85-CDI"; 

    } else if (!strcmp(keyword, "END")) {
      break;
    } else {
      sprintf(messageStr, "Unknown keyword: %s", keyword);
      Error(ERROR_LEVEL_1, routineName, messageStr);
    }

    if (nitem != Nitem) {
      sprintf(messageStr, "\n Read %d, not %d items (keyword = %s)\n",
	      nitem, Nitem, keyword);
      Error(ERROR_LEVEL_2, routineName, messageStr);
    }
    /* --- Spline interpolation in temperature T for all spatial
           locations. Linear if only 2 interpolation points given - - */

    //    if (strcmp(keyword, "TEMP") != 0) {
    if (!strcmp(keyword, "OMEGA") || !strcmp(keyword, "CE") ||
	!strcmp(keyword, "CI")    || !strcmp(keyword, "CP") ||
	!strcmp(keyword, "CH0")   || !strcmp(keyword, "CH+")||
	!strcmp(keyword, "CH") ) {

      if (Nitem > 2) {
	splineCoef(Nitem, T, coeff);
	splineEval(Nspace, atmos.T, C, hunt=TRUE);
      } else
	Linear(Nitem, T, coeff, Nspace, atmos.T, C, hunt=TRUE);
    }

    if (!strcmp(keyword, "OMEGA")) {

      /* --- Collisional excitation of ions --         -------------- */ 

      for (k = 0;  k < Nspace;  k++) {
        Cdown = C0 * atmos.ne[k] * C[k] /
                                 (atom->g[j] * sqrt(atmos.T[k]));
	atom->C[ij][k] += Cdown;
	atom->C[ji][k] += Cdown * atom->nstar[j][k]/atom->nstar[i][k];
      }
    } else if (!strcmp(keyword, "CE")) {      

      /* --- Collisional excitation of neutrals --     -------------- */ 

      gij = atom->g[i] / atom->g[j];
      for (k = 0;  k < Nspace;  k++) {
        Cdown = C[k] * atmos.ne[k] * gij * sqrt(atmos.T[k]);
	atom->C[ij][k] += Cdown;
	atom->C[ji][k] += Cdown * atom->nstar[j][k]/atom->nstar[i][k];
      }
    } else if (!strcmp(keyword, "CI")) {      

      /* --- Collisional ionization --                 -------------- */

      dE = atom->E[j] - atom->E[i];
      for (k = 0;  k < Nspace;  k++) {
        Cup = C[k] * atmos.ne[k] *
	  exp(-dE/(KBOLTZMANN*atmos.T[k])) * sqrt(atmos.T[k]);
	atom->C[ji][k] += Cup;
	atom->C[ij][k] += Cup * atom->nstar[i][k]/atom->nstar[j][k];
      }
    } else if (!strcmp(keyword, "CP")) {

      /* --- Collisions with protons --                -------------- */

      np = atmos.H->n[atmos.H->Nlevel-1];
      for (k = 0;  k < Nspace;  k++) {
        Cdown = np[k] * C[k];
	atom->C[ij][k] += Cdown;
	atom->C[ji][k] += Cdown * atom->nstar[j][k]/atom->nstar[i][k];
      }
    } else if (!strcmp(keyword, "CH")) {

      /* --- Collisions with neutral hydrogen --       -------------- */

      for (k = 0;  k < Nspace;  k++) {
        Cup = atmos.H->n[0][k] * C[k];
	atom->C[ji][k] += Cup;
	atom->C[ij][k] += Cup * atom->nstar[i][k]/atom->nstar[j][k];
      }
    } else if (!strcmp(keyword, "CH0")) {

      /* --- Charge exchange with neutral hydrogen --  -------------- */

      for (k = 0;  k < Nspace;  k++)
	atom->C[ij][k] += atmos.H->n[0][k] * C[k];

    } else if (!strcmp(keyword, "CH+")) {

      /* --- Charge exchange with protons --           -------------- */

      np = atmos.H->n[atmos.H->Nlevel-1];
      for (k = 0;  k < Nspace;  k++)
	atom->C[ji][k] += np[k] * C[k];

    }else if (!strcmp(keyword, "SHULL82")) {
      
      /* printf("shull82\n");*/
      
      acolsh=coeff[0];
      tcolsh=coeff[1];
      aradsh=coeff[2];
      xradsh=coeff[3];
      adish =coeff[4];
      bdish =coeff[5];
      t0sh  =coeff[6];
      t1sh  =coeff[7];
      
      for (k = 0;  k < Nspace;  k++) {

	summrs=sumscl*summers(i,j,atmos.ne[k],atom) + (1.0-sumscl) ;

	tg=atmos.T[k];
	
	cdn= aradsh*pow(tg/1.e4,-xradsh)  
	  + summrs * adish /tg/sqrt(tg) * exp(-t0sh/tg) * 
	  (1.0+bdish * (exp(- t1sh/tg)));
	
	cup=acolsh * sqrt(tg) * exp( -tcolsh / tg) / (1.0 + 0.1 * tg / tcolsh);
	
	/* convert coefficient from cm^3 s^-1 to m^3 s^-1 */
	cdn=cdn*pow(CM_TO_M,3);
	cup=cup*pow(CM_TO_M,3);

	cdn=cdn*atmos.ne[k];
	cup=cup*atmos.ne[k];

	/* 3-body recombination (high density limit) */
	cdn = cdn + cup * atom->nstar[i][k] / atom->nstar[j][k];
	
	atom->C[ij][k] += cdn;
	atom->C[ji][k] += cup;
      }
    }else if (!strcmp(keyword,"BADNELL")) {
      for (k = 0;  k < Nspace;  k++) {
	/*BADNELLL: Bhavna Rathore Jan-14 */
	/*fit for dielectronic recombination from Badnell
	  first line coefficients are the energies in K (ener in Chianti)
	  second line coefficients are the coefficients (coef in Chianti) */
       /*      printf("Badnell\n");*/
      	cdn=0.0;
	summrs=sumscl*summers(i,j,atmos.ne[k],atom) + (1.0-sumscl);
	tg=atmos.T[k];

        /* HU -- The outer loop shouldn't be here. It might work because
                 earlier a bigger, square, matrix was allocated, and the
                 additional elements that get addressed are
                 initialized to 0 (matrix uses calloc, rather than malloc)*/

	for (m = 0; m <= 1;  m++) {
          /* HU -- Nitem-1, or Nitem ? */

	  for (ii=0; ii< Nitem-1; ii++) {
	    cdn=cdn + badi[m+1][ii] * exp(-badi[m][ii] / tg);
	    
	    /* printf("\n m,n,badi[m][ii] %i,%i,%6.10lf",m,ii,badi[m][ii]);
	       printf("\n m+1,ii,badi[m+1][ii] %i,%i,%6.10lf",m+1,ii,badi[m+1][ii]); */
	  }
	}
	cdn=cdn * pow(tg,-1.5) ;
	/* convert coefficient from cm^3 s^-1 to m^3 s^-1 */
	cdn=cdn*pow(CM_TO_M,3);

        /* HU -- cup is not initialized, will be set to the values from some
	   previous entry so will create random values */

	cup=cup*pow(CM_TO_M,3);

	cdn=cdn* atmos.ne[k]*summrs;
	cup=cup*atmos.ne[k];

	/* 3-body recombination (high density limit) */
	cdn = cdn + cup * atom->nstar[i][k] / atom->nstar[j][k];
	
	atom->C[ij][k] += cdn;
	atom->C[ji][k] += cup;
      }
      /* HU -- matrix badi should be freed with freeMatrix((void **) badi) */

      /*badnell END    */ 
    }else if (!strcmp(keyword, "AR85-CDI")) {
      
      /*printf("ar85-cdi\n");*/
      
      /* Direct collionisional ionization */
      for (k = 0;  k < Nspace;  k++) {
	
	tg=atmos.T[k];
	cup=0.0;
	
	for (m = 0; m < Ncoef;  m++) {
	  
	  xj=cdi[m][0] * EV/KBOLTZMANN/tg;
	  fac=exp(-xj)*sqrt(xj);
	  
	  fxj= cdi[m][1]+cdi[m][2]*(1.0+xj) 
	    +(cdi[m][3]-xj*(cdi[m][1]+cdi[m][2] 
			    *(2.0+xj)))*fone(xj)+ 
	    cdi[m][4]*xj*ftwo(xj);
	  
	  fxj=fxj*fac;
	  fac = 6.69e-7 / pow(cdi[m][0],1.5);
	  cup = cup + fac*fxj * pow(CM_TO_M,3);
	  
	}
	
	if (cup<0) cup=0;

	cup=cup*atmos.ne[k];
	cdn=cup*atom->nstar[i][k]/atom->nstar[j][k];	  
	
	atom->C[ij][k] += cdn;
	atom->C[ji][k] += cup;
	
      }
      /* HU -- matrix cdi should be freed with freeMatrix((void **) cdi) */

    }else if (!strcmp(keyword,"AR85-CEA") ) {
          
      /* Autoionization */
      for (k = 0;  k < Nspace;  k++) {
	ar85cea(i,j,k,&fac, atom);
	cup=coeff[0]*fac*atmos.ne[k];
	atom->C[ji][k] += cup;
      }	  
      
    } else if (!strcmp(keyword, "AR85-CHP")) {
      
      /* charge transfer with ionized hydrogen */
      ar85t1=coeff[0];
      ar85t2=coeff[1];
      ar85a=coeff[2];
      ar85b=coeff[3];
      ar85c=coeff[4];
      ar85d=coeff[5];
      
      for (k = 0;  k < Nspace;  k++) {
	
	if (atmos.T[k]>=ar85t1 &&  atmos.T[k]<=ar85t2) {
	  
	  t4=atmos.T[k]/1.e4;
	  cup = ar85a * 1e-9 * pow(t4,ar85b) * exp(-ar85c*t4) 
	    * exp(-ar85d*EV/KBOLTZMANN/atmos.T[k])*atmos.H->n[5][k]
	    * pow(CM_TO_M,3);

	  atom->C[ji][k] += cup;
		  
	}
      }

  } else if (!strcmp(keyword, "AR85-CHH")) {
      
      /* charge transfer with neutral hydrogen */
      
      ar85t1=coeff[0];
      ar85t2=coeff[1];
      ar85a=coeff[2];
      ar85b=coeff[3];
      ar85c=coeff[4];
      ar85d=coeff[5];
      
      for (k = 0;  k < Nspace;  k++) {
	
	if (atmos.T[k]>=ar85t1 &&  atmos.T[k]<=ar85t2) {
	  
	  t4=atmos.T[k]/1.e4;
	  cdn = ar85a * 1e-9 * pow(t4,ar85b) * (1+ar85c*exp(ar85d*t4)) 
	    * atmos.H->n[0][k] * pow(CM_TO_M,3);

	  atom->C[ij][k] += cdn;
	  
	}
      }
      
    } else if (!strcmp(keyword, "BURGESS")) {
      
      /* Electron impact ionzation following Burgess & Chidichimo,
	 1982, MNRAS,203,1269-1280 */
      
      de= (atom->E[j]-atom->E[i]) / EV;
      zz=atom->stage[i];
      betab = 0.25 * ( sqrt( (100.0*zz +91.0) / (4.0*zz+3.0) ) -5.0 );
      cbar=2.3;
      
      for (k = 0;  k < Nspace;  k++) {
	
	
	dekt=de*EV/KBOLTZMANN/atmos.T[k];
	dekt=min(500,dekt);
	dekti=1.0/dekt;
        wlog=log(1.0 + dekti);
	wb= pow(wlog,betab/(1.0+dekti));
	cup = 2.1715e-8 * cbar * pow(13.6/de ,1.5) * sqrt(dekt) * E1(dekt) * wb *
	  atmos.ne[k] * pow(CM_TO_M,3);
	
	cup=cup*coeff[0]; //add fudge factor
	
	cdn = cup*atom->nstar[i][k]/atom->nstar[j][k];
	
	atom->C[ji][k] += cup;
	atom->C[ij][k] += cdn;
	
      }
      
    }
  }
  
  if (status == EOF) {
    sprintf(messageStr, "Reached end of datafile before all data was read");
    Error(ERROR_LEVEL_1, routineName, messageStr);
  }
  /* --- Clean up --                                   -------------- */

  free(C);
  free(T);
  free(coeff);

  sprintf(labelStr, "Collision Rate %2s", atom->ID);
  getCPU(3, TIME_POLL, labelStr);
}
/* ------- end ---------------------------- CollisionRate.c --------- */
