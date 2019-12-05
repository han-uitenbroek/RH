## ------- file: ------------------------- Makefile ----------------- ##
#
#      Version:       rh1.0, tools
#      Author:        Han Uitenbroek  (HUitenbroek@cfa.harvard.edu)
#      Last modified: Fri Oct 10 15:09:17 2003 --
#
##     --------------------------                      ----------RH-- ##

include ../makefiles/makefile.$(CPU).$(OS)


## --- If no FORTRAN compiler is available remove -lrh_f90 in following

    LIBS = -lrh -lrh_f90 $(ARCHLIBS)


## --- Rules for the executables --                    -------------- ##

all:     addgrad  convertatom  listincl  topbase  printneff  impact \
         hcollisions  make_h  listincl_f  avgmolweight  readnist  nist_pf \
         waveinfo  readkurucz


addgrad:  addgrad.o  modelatom.o  rawatom.o librh
	$(CC) -o $@  $(LDFLAGS) addgrad.o  modelatom.o rawatom.o \
         $(LIBS)

convertatom:  convertatom.o modelatom.o librh
	$(CC) -o $@  $(LDFLAGS) convertatom.o modelatom.o $(LIBS)

topbase:  topbase.o modelatom.o librh
	$(CC) -o $@  $(LDFLAGS) topbase.o modelatom.o $(LIBS)

#-# adjustBF:  adjustbf.o rawatom.o modelatom.o librh
#-#	$(CC) -o $@  $(LDFLAGS) adjustbf.o modelatom.o \
#-#         rawatom.o $(LIBS)

listincl:  listincl.o
	$(CC) -o $@  $(LDFLAGS) listincl.o

listincl_f:  listincl_f.o
	$(CC) -o $@  $(LDFLAGS) listincl_f.o

printneff:  printneff.o rawatom.o librh
	$(CC) -o $@  $(LDFLAGS) printneff.o rawatom.o $(LIBS)

impact:  impact.o rawatom.o laguerre.o bessel.o librh
	$(CC) -o $@  $(LDFLAGS) impact.o rawatom.o \
         laguerre.o bessel.o $(LIBS)

hcollisions:  hcollisions.o laguerre.o librh
	$(CC) -o $@  $(LDFLAGS) hcollisions.o laguerre.o $(LIBS)

make_h:  make_h.o modelatom.o librh
	$(CC) -o $@  $(LDFLAGS) make_h.o modelatom.o \
        $(LIBS)

avgmolweight: avgmolweight.o librh
	$(CC) -o $@  $(LDFLAGS) avgmolweight.o $(LIBS)

readnist: readnist.o librh
	$(CC) -o $@  $(LDFLAGS) readnist.o $(LIBS)

nist_pf: nist_pf.o librh
	$(CC) -o $@  $(LDFLAGS) nist_pf.o $(LIBS)

waveinfo: waveinfo.o librh
	$(CC) -o $@  $(LDFLAGS) waveinfo.o $(LIBS)

constructatom: constructatom.o modelatom.o librh
	$(CC) -o $@  $(LDFLAGS) constructatom.o modelatom.o $(LIBS)

readkurucz:  readkurucz.o writemulti.o
	$(CC) -o $@  $(LDFLAGS) readkurucz.o writemulti.o $(LIBS)

readtrujillo:  readtrujillo.o writemulti.o
	$(CC) -o $@  $(LDFLAGS) readtrujillo.o writemulti.o $(LIBS)

##lande:  lande.o rawatom.o librh
##	$(CC) -o $@  $(LDFLAGS) lande.o rawatom.o $(LIBS)


## --- If no FORTRAN compiler is available remove librh_f90.a in following

librh:
	cd ..; $(MAKE) librh.a librh_f90.a


## --- Clean up --                                     -------------- ##

clean:
	rm -f  *.o  core  addgrad  convertatom  listincl  topbase  printneff \
          make_h  hcollisions  impact  listincl_f  avgmolweight  \
          readnist  nist_pf  waveinfo  readkurucz  readtrujillo


## --- Explicit dependencies on include files --       -------------- ##


addgrad.o:              ../rh.h  ../atom.h  ../error.h  \
                        ../inputs.h

adjustbf.o:             ../rh.h  ../atom.h  ../constant.h  \
                        ../error.h  ../inputs.h

avgmolweight.o:         ../rh.h  ../atom.h  ../atmos.h  \
                        ../inputs.h

bessel.o:               
constructatom.o:        ../rh.h  ../atom.h  ../atmos.h  \
                        ../spectrum.h  ../constant.h  ../inputs.h \
                        ../error.h  \
                        ../statistics.h  ../atomweights.h

convertatom.o:          ../rh.h  ../atom.h  ../constant.h  \
                        ../error.h  ../inputs.h

hcollisions.o:          ../rh.h  ../constant.h  ../error.h  \
                        ../inputs.h

impact.o:               ../rh.h  ../atom.h  ../constant.h  \
                        ../error.h  ../inputs.h

laguerre.o:             
lande.o:                ../rh.h  ../atom.h  ../error.h  \
                        ../inputs.h

listincl.o:             
listincl_f.o:           
make_h.o:               ../rh.h  ../atom.h  ../constant.h  \
                        ../error.h  ../inputs.h

modelatom.o:            ../rh.h  ../atom.h  ../constant.h  \
                        ../error.h  ../inputs.h

nist_pf.o:              ../rh.h  ../atom.h  ../constant.h  \
                        ../inputs.h  ../error.h

printneff.o:            ../rh.h  ../atom.h  ../constant.h  \
                        ../error.h  ../inputs.h

rawatom.o:              ../rh.h  ../atom.h  ../constant.h  \
                        ../error.h  ../inputs.h  ../atomweights.h

readnist.o:             ../rh.h  ../atom.h  ../constant.h  \
                        ../inputs.h  ../error.h

readkurucz.o:           ../rh.h  ../atom.h  ../atomweights.h  ../atmos.h \
                        ../rhf1d/geometry.h  ../error.h  ../inputs.h \
                        ../constant.h

readtrujillo.o:         ../rh.h  ../atom.h  ../atmos.h  ../rhf1d/geometry.h \
                        ../error.h  ../inputs.h  ../constant.h  

topbase.o:              ../rh.h  ../atom.h  ../constant.h  \
                        ../error.h  ../inputs.h

waveinfo.o:             

writemulti.o:           ../rh.h  ../atom.h  ../atmos.h  ../rhf1d/geometry.h \
                        ../error.h  ../constant.h  ../xdr.h  

## ------- end ---------------------------- Makefile ---------------- ##
