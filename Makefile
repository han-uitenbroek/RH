## ------- file: ------------------------- Makefile ----------------- ##
#
#      Version:       rh2.0
#      Author:        Han Uitenbroek (huitenbroek@nso.edu)
#      Last modified: Thu May 24 16:09:35 2018 --
#
##     --------------------------                      ----------RH-- ##

include makefile.$(CPU).$(OS)


all: librh.a(getcpu.o) librh.a(fpehandler.o) librh.a


## --- Rules for the library --                        -------------- ##

librh.a: \
 librh.a(abundance.o) \
 librh.a(accelerate.o) \
 librh.a(background.o) \
 librh.a(backgropac_xdr.o) \
 librh.a(barklem.o) \
 librh.a(bezier_aux.o) \
 librh.a(broad.o) \
 librh.a(brs_xdr.o) \
 librh.a(chemequil.o) \
 librh.a(cocollisions.o) \
 librh.a(collision.o) \
 librh.a(complex.o) \
 librh.a(cubeconvol.o) \
 librh.a(duplicate.o) \
 librh.a(error.o) \
 librh.a(expint.o) \
 librh.a(expspline.o) \
 librh.a(fillgamma.o) \
 librh.a(fixedrate.o) \
 librh.a(fpehandler.o) \
 librh.a(gammafunc.o) \
 librh.a(gaussleg.o) \
 librh.a(getcpu.o) \
 librh.a(getlambda.o) \
 librh.a(getline.o) \
 librh.a(giigen.o) \
 librh.a(h2collisions.o) \
 librh.a(hunt.o) \
 librh.a(humlicek.o) \
 librh.a(hydrogen.o) \
 librh.a(initial_xdr.o) \
 librh.a(initscatter.o) \
 librh.a(iterate.o) \
 librh.a(kurucz.o) \
 librh.a(linear.o) \
 librh.a(ltepops.o) \
 librh.a(ludcmp.o) \
 librh.a(matrix.o) \
 librh.a(maxchange.o) \
 librh.a(metal.o) \
 librh.a(molzeeman.o) \
 librh.a(nemetals.o) \
 librh.a(ohchbf.o) \
 librh.a(opacity.o) \
 librh.a(options.o) \
 librh.a(order.o) \
 librh.a(parse.o) \
 librh.a(paschen.o) \
 librh.a(planck.o) \
 librh.a(pops_xdr.o) \
 librh.a(profile.o) \
 librh.a(radrate_xdr.o) \
 librh.a(readatom.o) \
 librh.a(readb_xdr.o) \
 librh.a(readj.o) \
 librh.a(readvalue.o) \
 librh.a(rayleigh.o) \
 librh.a(readinput.o) \
 librh.a(readmolecule.o) \
 librh.a(redistribute.o) \
 librh.a(scatter.o) \
 librh.a(solvene.o) \
 librh.a(sortlambda.o) \
 librh.a(spline.o) \
 librh.a(statequil.o) \
 librh.a(stokesopac.o) \
 librh.a(stopreq.o) \
 librh.a(thomson.o) \
 librh.a(vacuumtoair.o) \
 librh.a(voigt.o) \
 librh.a(w3.o) \
 librh.a(wigner.o) \
 librh.a(writeatmos_xdr.o) \
 librh.a(writeatom_xdr.o) \
 librh.a(writecoll_xdr.o) \
 librh.a(writedamp_xdr.o) \
 librh.a(writeinput_xdr.o) \
 librh.a(writemetal_xdr.o) \
 librh.a(writemolec_xdr.o) \
 librh.a(writeopac_xdr.o) \
 librh.a(writespect_xdr.o) \
 librh.a(zeeman.o)


## FORTRAN 90 alternatives --                          -------------- ##

librh_f90.a: \
 librh_f90.a(hui_.o) \
 librh_f90.a(humlicek_.o)


## --- Specific compilation rules for machine-dependent files -- ---- ##

librh.a(fpehandler.o):  fpehandler.c
	$(CC) -c $(CFLAGS) -D$(CPU) -D$(OS)  $*.c
	ar $(ARFLAGS) librh.a $%
	rm -f $%

librh.a(getcpu.o):  getcpu.c
	$(CC) -c $(CFLAGS) -D$(CPU) -D$(OS)  $*.c
	ar $(ARFLAGS) librh.a $%
	rm -f $%


## --- Clean up --                                     -------------- ##

clean:
	rm -f *.o  *.a


## --- Explicit dependencies on include files --       -------------- ##


librh.a(abundance.o):     rh.h  atom.h  atmos.h  constant.h  error.h  inputs.h  xdr.h  atomweights.h
librh.a(accelerate.o):    rh.h  accelerate.h  error.h  statistics.h
librh.a(backgropac_xdr.o):rh.h  atom.h  atmos.h  spectrum.h  background.h  error.h  inputs.h  xdr.h
librh.a(background.o):    rh.h  atom.h  atmos.h  spectrum.h  constant.h  background.h  error.h  statistics.h  inputs.h  xdr.h
librh.a(barklem.o):       rh.h  atom.h  atmos.h  constant.h  error.h
librh.a(bezier_aux.o):    bezier.h
librh.a(broad.o):         rh.h  atom.h  atmos.h  constant.h  error.h
librh.a(brs_xdr.o):       rh.h  atom.h  atmos.h  spectrum.h  error.h  xdr.h
librh.a(chemequil.o):     rh.h  atom.h  atmos.h  background.h  accelerate.h  constant.h  error.h  statistics.h  inputs.h
librh.a(cocollisions.o):  rh.h  atom.h  atmos.h  error.h  constant.h
librh.a(collision.o):     rh.h  atom.h  atmos.h  constant.h  error.h  inputs.h  statistics.h
librh.a(complex.o):       rh.h  complex.h
librh.a(cubeconvol.o):    rh.h
librh.a(duplicate.o):     rh.h  atom.h  background.h  error.h
librh.a(error.o):         rh.h  error.h  inputs.h
librh.a(expint.o):        rh.h  error.h
librh.a(expspline.o):     rh.h
librh.a(fillgamma.o):     rh.h  error.h  atom.h  atmos.h  spectrum.h  atom.h  inputs.h  constant.h
librh.a(fixedrate.o):     rh.h  atom.h  atmos.h  constant.h  statistics.h
librh.a(fpehandler.o):    rh.h  error.h
librh.a(gammafunc.o):     rh.h
librh.a(gaussleg.o):      rh.h  constant.h
librh.a(getcpu.o):        rh.h  statistics.h  error.h  inputs.h
librh.a(getlambda.o):     rh.h  atom.h  atmos.h  constant.h  error.h  inputs.h
librh.a(getline.o):       rh.h  error.h
librh.a(giigen.o):        rh.h  atom.h  constant.h
librh.a(cocollisions.o):  rh.h  atom.h  atmos.h  error.h  constant.h
librh.a(humlicek.o):      rh.h  complex.h
librh.a(hunt.o):          rh.h
librh.a(hydrogen.o):      rh.h  atom.h  atmos.h  constant.h  background.h  error.h
librh.a(initial_xdr.o):   rh.h  atom.h  atmos.h  spectrum.h  accelerate.h  constant.h  statistics.h  error.h  inputs.h  xdr.h
librh.a(initscatter.o):   rh.h  atom.h  atmos.h  spectrum.h  accelerate.h  error.h  statistics.h  inputs.h
librh.a(iterate.o):       rh.h  atom.h  atmos.h  spectrum.h  background.h  accelerate.h  error.h  statistics.h  inputs.h
librh.a(kurucz.o):        rh.h  atom.h  atmos.h  background.h  spectrum.h  constant.h  inputs.h  error.h
librh.a(linear.o):        rh.h
librh.a(ltepops.o):       rh.h  atom.h  atmos.h  constant.h  error.h  statistics.h
librh.a(ludcmp.o):        rh.h  error.h
librh.a(matrix.o):        rh.h  error.h
librh.a(maxchange.o):     rh.h  accelerate.h  error.h  inputs.h
librh.a(metal.o):         rh.h  atom.h  atmos.h  error.h  constant.h  background.h  inputs.h
librh.a(molzeeman.o):     rh.h  atom.h  spectrum.h  error.h
librh.a(nemetals.o):      rh.h  atom.h  atmos.h  background.h
librh.a(ohchbf.o):        rh.h  atom.h  atmos.h  constant.h
librh.a(opacity.o):       rh.h  error.h  atom.h  atmos.h  spectrum.h  inputs.h  constant.h
librh.a(options.o):       rh.h  error.h  inputs.h
librh.a(order.o):         rh.h
librh.a(parse.o):         rh.h  error.h  inputs.h
librh.a(paschen.o):       rh.h  error.h
librh.a(planck.o):        rh.h  atom.h  spectrum.h  constant.h
librh.a(pops_xdr.o):      rh.h  atom.h  atmos.h  error.h  xdr.h
librh.a(profile.o):       rh.h  atom.h  atmos.h  inputs.h  constant.h  statistics.h  error.h
librh.a(radrate_xdr.o):   rh.h  atom.h  atmos.h  error.h  inputs.h  xdr.h
librh.a(rayleigh.o):      rh.h  atom.h  atmos.h  constant.h  background.h  error.h
librh.a(readatom.o):      rh.h  atom.h  atmos.h  spectrum.h  background.h  constant.h  error.h  inputs.h  statistics.h
librh.a(readb_xdr.o):     rh.h  atom.h  atmos.h  error.h  inputs.h  xdr.h
librh.a(readinput.o):     rh.h  atom.h  atmos.h  spectrum.h  constant.h  statistics.h  inputs.h  error.h
librh.a(readj.o):         rh.h  atom.h  atmos.h  spectrum.h  background.h  inputs.h  error.h
librh.a(readmolecule.o):  rh.h  atom.h  atmos.h  spectrum.h  background.h  constant.h  error.h  inputs.h  statistics.h
librh.a(readvalue.o):     rh.h  atom.h  atmos.h  error.h  inputs.h
librh.a(redistribute.o):  rh.h  atom.h  atmos.h  spectrum.h  accelerate.h  error.h  inputs.h
librh.a(scatter.o):       rh.h  atom.h  atmos.h  spectrum.h  inputs.h  constant.h  error.h  statistics.h
librh.a(solvene.o):       rh.h  atom.h  atmos.h  constant.h  error.h  statistics.h
librh.a(sortlambda.o):    rh.h  atom.h  atmos.h  spectrum.h  background.h  constant.h  inputs.h  error.h  statistics.h  xdr.h
librh.a(spline.o):        rh.h
librh.a(statequil.o):     rh.h  atom.h  atmos.h  inputs.h  accelerate.h  statistics.h  error.h
librh.a(stokesopac.o):    rh.h  atom.h  atmos.h  spectrum.h  inputs.h
librh.a(stopreq.o):       rh.h  error.h
librh.a(thomson.o):       rh.h  atom.h  atmos.h  constant.h  background.h
librh.a(vacuumtoair.o):   rh.h  atom.h  spectrum.h
librh.a(voigt.o):         rh.h  constant.h  complex.h  error.h
librh.a(w3.o):            rh.h
librh.a(wigner.o):        rh.h
librh.a(writeatmos_xdr.o):rh.h  atom.h  atmos.h  error.h  inputs.h  xdr.h
librh.a(writeatom_xdr.o): rh.h  atom.h  atmos.h  spectrum.h  inputs.h  constant.h  error.h  xdr.h
librh.a(writecoll_xdr.o): rh.h  atom.h  atmos.h  error.h  inputs.h  xdr.h
librh.a(writedamp_xdr.o): rh.h  atom.h  atmos.h  error.h  inputs.h  xdr.h
librh.a(writeinput_xdr.o):rh.h  atom.h  atmos.h  spectrum.h  error.h  inputs.h  xdr.h
librh.a(writemetal_xdr.o):rh.h  atom.h  atmos.h  error.h  xdr.h
librh.a(writemolec_xdr.o):rh.h  atom.h  atmos.h  error.h  xdr.h
librh.a(writeopac_xdr.o): rh.h  atom.h  atmos.h  spectrum.h  error.h  inputs.h  xdr.h
librh.a(writespect_xdr.o):rh.h  atom.h  atmos.h  spectrum.h  constant.h  error.h  inputs.h  xdr.h
librh.a(zeeman.o):        rh.h  atom.h  atmos.h  error.h  inputs.h  statistics.h

## ------- end ---------------------------- Makefile ---------------- ##
