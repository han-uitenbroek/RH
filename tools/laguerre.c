/* ------- file: -------------------------- laguerre.c --------------

       Version:       rh1.0, tools
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Oct 8, 1996

       --------------------------                      ----------RH-- */

/* --- Routine for integration with Gauss-Laguerre for integrals of type:

       /oo
       | f(x) exp(-x) dx
       /0

       --                                              -------------- */

 
#define N_LAGUERRE 8

/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */


/* ------- begin -------------------------- GaussLaguerre.c --------- */

double GaussLaguerre( double (*function) (double x) )
{
  register int nl;

  double integral = 0.0;

  static double x_Laguerre[N_LAGUERRE] =
    { 0.170279632305,  0.903701776799,  2.251086629866,  4.266700170288,
      7.045905402393, 10.758516010181, 15.740678641278, 22.863131736889 };
  static double w_Laguerre[N_LAGUERRE] =
    { 3.69188589342E-01, 4.18786780814E-01, 1.75794986637E-01,
      3.33434922612E-02, 2.79453623523E-03, 9.07650877336E-05,
      8.48574671627E-07, 1.04800117487E-09 };

  for (nl = 0;  nl < N_LAGUERRE;  nl++)
    integral += w_Laguerre[nl] * function(x_Laguerre[nl]);

  return integral;
}
/* ------- end ---------------------------- GaussLaguerre.c --------- */
