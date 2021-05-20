: ------------------------- Makefile ----------------- ##
#
#      Version:       rh2.0, 2-D Cartesian
#      Author:        Han Uitenbroek (huitenbroek@nso.edu)
#      Last modified: Thu May 20 13:06:17 2021 --
#
##     --------------------------                      ----------RH-- ##

GEOMETRY = TWO_D_PLANE

include ../makefiles/rhmakefile


## --- Define groups of object files --                -------------- ##


TWO_D_OBJS = anglequad.o  bezier_2D.o  formal.o  geometry.o \
             readatmos_xdr.o  rhsc2d.o  piecewise_2D.o  piecestokes_2D.o \
             solvelong.o  project.o  rii.o  writegeom_xdr.o  writeflux_xdr.o

   TWODOBJ = anglequad.o  bezier_2D.o  iter_2d.o  geometry.o \
             piecewise_2D.o  piecestokes_2D.o  solve2d.o  \
             solvelong.o  project.o

  FORMOBJS = anglequad.o  formal2d.o  geometry.o \
             piecewise_2D.o  piecestokes_2D.o  solvelong.o  project.o

  RAY_OBJS = anglequad.o  bezier_2D.o  formal.o  geometry.o  rii.o \
             readatmos_xdr.o  solveray.o  piecewise_2D.o  piecestokes_2D.o \
             solvelong.o  project.o

CHART_OBJS = anglequad.o  bezier_2D.o  chartflow.o  formal.o  geometry.o \
             readatmos_xdr.o  rii.o  piecewise_2D.o  piecestokes_2D.o \
             solvelong.o  project.o

BACKGROPAC = anglequad.o  backgrcontr.o  readatmos_xdr.o  project.o


## --- Rules for the executables --                    -------------- ##

all:   rhsc2d  solve2d  formal2d  solveray  chartflow  backgrcontr

rhsc2d:  $(TWO_D_OBJS)  librh
	$(LD) -o $@  $(LDFLAGS) $(TWO_D_OBJS) $(LIBS)

solve2d:  $(TWODOBJ)  librh
	$(LD) -o $@  $(LDFLAGS) $(TWODOBJ) $(LIBS)

formal2d:  $(FORMOBJS)  librh
	$(LD) -o $@  $(LDFLAGS) $(FORMOBJS)  $(LIBS)

solveray:  $(RAY_OBJS)  librh
	$(LD) -o $@  $(LDFLAGS) $(RAY_OBJS) $(LIBS)

chartflow:  $(CHART_OBJS)  librh
	$(LD) -o $@  $(LDFLAGS) $(CHART_OBJS) $(LIBS)

backgrcontr:  $(BACKGROPAC) librh
	$(LD) -o $@  $(LDFLAGS) $(BACKGROPAC) $(LIBS)


## --- If no FORTRAN compiler is available remove librh_f90.a in following

librh:
	cd ..; $(MAKE) librh.a librh_f90.a


## --- Clean up --                                     -------------- ##

clean:
	rm -f *.o  rhsc2d  solve2d  formal2d  solveray  chartflow  backgrcontr


## --- Explicit dependencies on include files --       -------------- ##

anglequad.o:            ../rh.h  ../atom.h  ../atmos.h  \
                        geometry.h  ../inputs.h  ../error.h  ../constant.h

bezier_2D.o:            ../rh.h  ../atom.h  ../atmos.h       \
                        geometry.h  ../spectrum.h  ../bezier.h

backgrcontr.o:          ../rh.h  ../atom.h  ../atmos.h  \
                        ../spectrum.h  geometry.h  ../background.h  \
                        ../inputs.h  ../error.h  ../statistics.h

chartflow.o:            ../rh.h  ../atom.h  ../atmos.h  \
                        geometry.h  ../spectrum.h  ../background.h  \
                        ../statistics.h  ../inputs.h  ../error.h  ../xdr.h

formal.o:               ../rh.h  ../atom.h  ../atmos.h  \
                        geometry.h  ../spectrum.h  ../error.h  \
                        ../background.h  ../inputs.h  ../constant.h  \
                        ../statistics.h  ../xdr.h

formal2d.o:             ../rh.h  ../atom.h  ../atmos.h  ../spectrum.h \
                        geometry.h  ../inputs.h  ../error.h  \
                        ../statistics.h  ../constant.h

geometry.o:             ../rh.h  ../atom.h  ../atmos.h  \
                        geometry.h  ../inputs.h  ../statistics.h  ../error.h

iter_2d.o:              ../rh.h  ../atom.h  ../atmos.h  \
                        geometry.h  ../error.h  ../accelerate.h  \
                        ../statistics.h

project.o:              ../rh.h  ../atom.h  ../atmos.h  \
                        geometry.h

readatmos_xdr.o:        ../rh.h  ../atom.h  ../atmos.h  \
                        geometry.h  ../spectrum.h  ../background.h  \
                        ../constant.h  ../statistics.h  ../error.h  \
                        ../inputs.h  ../xdr.h

rhsc2d.o:               ../rh.h  ../atom.h  ../atmos.h  \
                        geometry.h  ../spectrum.h  ../background.h  \
                        ../statistics.h  ../inputs.h  ../xdr.h

rii.o:                  ../rh.h  ../atom.h  ../atmos.h  \
                        geometry.h  ../constant.h  ../error.h

piecewise_2D.o:         ../rh.h  ../atom.h  ../atmos.h  ../spectrum.h \
                        geometry.h  ../inputs.h  ../error.h

piecestokes_2D.o:       ../rh.h  ../atom.h  ../atmos.h  ../spectrum.h \
                        geometry.h  ../inputs.h  ../error.h

solve15d.o:             ../rh.h  ../atom.h  ../atmos.h  \
                        geometry.h  ../spectrum.h  ../background.h  \
                        ../accelerate.h  ../error.h  ../constant.h  \
                        ../inputs.h  ../statistics.h

solve2d.o:              ../rh.h  ../atom.h  ../atmos.h  ../spectrum.h \
                        geometry.h  ../accelerate.h  ../inputs.h  \
                        ../error.h  ../statistics.h  ../constant.h  \
                        ../xdr.h

solvelong.o:            ../rh.h  ../atom.h  ../atmos.h  ../spectrum.h \
                        geometry.h

solveray.o:             ../rh.h  ../atom.h  ../atmos.h  \
                        geometry.h  ../spectrum.h  ../background.h  \
                        ../statistics.h  ../inputs.h  ../error.h  ../xdr.h

writeflux_xdr.o:        ../rh.h  ../atom.h  ../atmos.h  \
                        geometry.h  ../spectrum.h  ../constant.h  \
                        ../error.h  ../xdr.h

writegeom_xdr.o:        ../rh.h  ../atom.h  ../atmos.h  \
                        geometry.h  ../error.h  ../inputs.h  ../xdr.h

## ------- end ---------------------------- Makefile ---------------- ##
