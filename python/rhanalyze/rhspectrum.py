import sys

if sys.version_info[1] < 13:
    import xdrlib
else:
    import mda_xdrlib.xdrlib as xdrlib

import numpy as np

import rhanalyze.rhgeometry
from rhanalyze.rhtools import read_farray, read_string
from rhanalyze.rhtools import vacuum_to_air, air_to_vacuum

class spectrum:
    
    def __init__(self, inputs, geometry, atmos, filename='spectrum.out', \
                 fluxfile='flux.out'):

        self.filename = filename
        self.fluxfile = fluxfile
        self.read(inputs, geometry, atmos)
        self.read_flux(geometry)


    def read(self, inputs, geometry, atmos):

        f  = open(self.filename, 'rb')
        up = xdrlib.Unpacker(f.read())
        f.close()

        self.Nspect = up.unpack_int()
        self.waves = read_farray(self.Nspect, up, "double")

        if geometry.type == "ONE_D_PLANE" or\
           geometry.type == 'SPHERICAL_SYMMETRIC':
            dim = [geometry.Nrays, self.Nspect]

        elif geometry.type == 'TWO_D_PLANE':
            dim = [geometry.Nx, geometry.Nrays, self.Nspect]
        elif geometry.type == 'THREE_D_PLANE':
            dim = [geometry.Nx, geometry.Ny, geometry.Nrays, self.Nspect]

        self.I = read_farray(dim, up, "double")

        self.vacuum_to_air = up.unpack_int()
        self.air_limit     = float(up.unpack_double())

        if atmos.stokes or inputs.backgr_pol:
            self.Q = read_farray(dim, up, "double")
            self.U = read_farray(dim, up, "double")
            self.V = read_farray(dim, up, "double")

        up.done()


    def read_flux(self, geometry):

        f  = open(self.fluxfile, 'rb')
        up = xdrlib.Unpacker(f.read())
        f.close()
        
        match geometry.type:
            case "ONE_D_PLANE" | "SPHERICAL_SYMMETRIC":
                dim = [self.Nspect]
            case "TWO_D_PLANE":
                dim = [geometry.Nx, self.Nspect]
            case "THREE_D_PLANE":
                dim = [geometry.Nx, geometry.Ny, self.Nspect]

        self.flux = read_farray(dim, up, "double")
        up.done()


class rays:

    def __init__(self, inputs, geometry, spectrum, filename='spectrum_1.00'):
        self.filename = filename
        self.read(inputs, geometry, spectrum)
 
    def read(self, inputs, geometry, spectrum):

        f  = open(self.filename, 'rb')
        up = xdrlib.Unpacker(f.read())
        f.close()
        
        if geometry.type == "ONE_D_PLANE":
            self.muz = up.unpack_double()
            
            dim1 = [spectrum.Nspect]
            dim2 = [geometry.Ndep]
            
        elif geometry.type == 'SPHERICAL_SYMMETRIC':
            self.muz = up.unpack_double()
           
            dim1 = [spectrum.Nspect]
            dim2 = [geometry.Nradius]
            
        elif geometry.type == 'TWO_D_PLANE':
            self.mux = up.unpack_double()
            self.muz = up.unpack_double()

            dim1 = [geometry.Nx, spectrum.Nspect]
            dim2 = [geometry.Nx, geometry.Nz]

        elif geometry.type == 'THREE_D_PLANE':
            self.mux = up.unpack_double()
            self.muy = up.unpack_double()

            dim1 = [geometry.Nx, geometry.Ny, spectrum.Nspect]
            dim2 = [geometry.Nx, geometry.Ny, geometry.Nz]


        self.I = read_farray(dim1, up, "double")

        self.Nopac = up.unpack_int()
        if self.Nopac > 0:
            self.opac = {}
            for n in range(self.Nopac):
                self.opac[n] = opac(dim2, up)

        try: spectrum.Q
        except AttributeError:
            up.done()
            return
        else:
            self.Q = read_farray(dim1, up, "double")
            self.U = read_farray(dim1, up, "double")
            self.V = read_farray(dim1, up, "double")

            up.done()


class opac:

    def __init__(self, dim, up):
        self.read(dim, up)
        
    def read(self, dim, up):
        self.nspect = up.unpack_int()
        self.chi    = read_farray(dim, up, "double")
        self.S      = read_farray(dim, up, "double")
    
