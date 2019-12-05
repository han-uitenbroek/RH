## ------- file: ------------------------- Makefile ----------------- ##
#
#      Version:       rh2.0, 1-D spherically symmetric
#      Author:        Han Uitenbroek (huitenbroek@nso.edu)
#      Last modified: Wed Apr 22 09:41:23 2009 --
#
##     --------------------------                      ----------RH-- ##


GEOMETRY = SPHERICAL_SYMMETRIC

include ../makefiles/makefile.$(CPU).$(OS)


## --- If no FORTRAN compiler is available remove -lrh_f90 in following

    LIBS = -lrh -lrh_f90 $(ARCHLIBS)


## --- Define groups of object files --                -------------- ##

SPHEREOBJS = feautrier.o  formal.o  geometry.o  readatmos.o  rhsphere.o \
             rii.o  project.o  writegeom_xdr.o


## --- Rules for the executables --                    -------------- ##


rhsphere:  $(SPHEREOBJS)  librh
	$(LD) -o $@  $(LDFLAGS) $(SPHEREOBJS) $(LIBS)


## --- If no FORTRAN compiler is available remove librh_f90.a in following

librh:
	cd ..; $(MAKE) librh.a  librh_f90.a


## --- Clean up --                                     -------------- ##

clean:
	rm -f *.o  rhsphere


## --- Explicit dependencies on include files --       -------------- ##


feautrier.o:            ../rh.h  ../error.h  ../atom.h  ../atmos.h  \
                        geometry.h  ../spectrum.h  ../error.h

formal.o:               ../rh.h  ../atom.h  ../atmos.h  \
                        ../spectrum.h  geometry.h  ../error.h  ../inputs.h  \
                        ../background.h  ../constant.h  ../statistics.h  \
                        ../xdr.h

geometry.o:             ../rh.h  ../atom.h  ../atmos.h  \
                        geometry.h  ../statistics.h

project.o:              ../rh.h  ../atom.h  ../atmos.h  \
                        geometry.h

readatmos.o:            ../rh.h  ../atom.h  ../atmos.h  \
                        geometry.h  ../spectrum.h  ../background.h  \
                        ../constant.h  ../statistics.h  ../error.h  \
                        ../inputs.h  ../xdr.h

readinput.o:            ../rh.h  ../atmos.h  ../spectrum.h  \
                        ../constant.h  ../statistics.h  ../inputs.h  \
                        ../error.h

rhsphere.o:             ../rh.h  ../atom.h  ../atmos.h  \
                        geometry.h  ../spectrum.h  ../background.h  \
                        ../inputs.h  ../statistics.h  ../xdr.h

rii.o:                  ../rh.h  ../atom.h  ../error.h

solveplane.o:           ../rh.h  ../atom.h  ../atmos.h  \
                        geometry.h  ../spectrum.h  ../background.h  \
                        ../accelerate.h  ../constant.h  ../inputs.h \
                        ../error.h  ../statistics.h

writegeom_xdr.o:        ../rh.h  ../atom.h  ../atmos.h  \
                        geometry.h  ../error.h  ../inputs.h  ../xdr.h

## ------- end ---------------------------- Makefile ---------------- ##
