import xdrlib
import numpy as np
import os.path

import rhanalyze.rhgeometry
from rhanalyze.rhtools import read_farray, read_string

class atoms:

    def __init__(self, geometry, path='./atom.H.out'):
        self.filename = path
        self.read(geometry)

        (dirname, filename) = os.path.split(path)
        
        popsfile = '{0}/pops.{1}.out'.format(dirname, self.atomID)
        if os.path.exists(popsfile):
            self.readpops(geometry, popsfile)

        ratesfile = '{0}/radrate.{1}.out'.format(dirname, self.atomID)
        if os.path.exists(ratesfile):
            self.readrates(geometry, ratesfile)

        dampingfile = '{0}/damping.{1}.out'.format(dirname, self.atomID)
        if os.path.exists(dampingfile):
            self.readdamping(geometry, dampingfile)

        collisionfile = '{0}/collrate.{1}.out'.format(dirname, self.atomID)
        if os.path.exists(collisionfile):
            self.readcollisions(geometry, collisionfile)
            
    def read(self, geometry):

        f  = open(self.filename, 'rb')
        up = xdrlib.Unpacker(f.read())
        f.close()

        self.active = up.unpack_int()
        self.Nlevel = up.unpack_int()
        self.Nline  = up.unpack_int()
        self.Ncont  = up.unpack_int()
        self.Nfixed = up.unpack_int()

        self.abund  = up.unpack_double()
        self.weight = up.unpack_double()

        self.labels = {}
        for i in range(self.Nlevel):
            self.labels[i] = read_string(up)
            
        self.atomID = self.labels[0][0:2].strip()

        self.g     = read_farray([self.Nlevel], up, "double")
        self.E     = read_farray([self.Nlevel], up, "double")
        self.stage = read_farray([self.Nlevel], up, "int")

        Nrad = self.Nline + self.Ncont

        self.transition = {}
        for kr in range(Nrad):
            self.transition[kr] = atoms.transition(up)

        for kr in range(self.Nline, Nrad):
            if self.transition[kr].shape == "HYDROGENIC":
                self.transition[kr].waves = up.unpack_double()
            else:
                self.transition[kr].waves =\
                read_farray([self.transition[kr].Nwave], up, "double")
                self.transition[kr].alpha =\
                read_farray([self.transition[kr].Nwave], up, "double")
            
        self.fixed = {}
        for kr in range(self.Nfixed):
            self.transition[kr] = atoms.fixed(up)

        up.done()

 
    def readpops(self, geometry, popsfile):
        
        f  = open(popsfile, 'rb')
        up = xdrlib.Unpacker(f.read())
        f.close()

        atmosID     = read_string(up)
        Nlevel      = up.unpack_int()
        Nspace      = up.unpack_int()

        if geometry.type == "ONE_D_PLANE":
            dim = [geometry.Ndep, self.Nlevel]
        elif geometry.type == 'SPHERICAL_SYMMETRIC':
            dim = [geometry.Nradius, self.Nlevel]
        elif geometry.type == 'TWO_D_PLANE':
            dim = [geometry.Nx, geometry.Nz, self.Nlevel]
        elif geometry.type == 'THREE_D_PLANE':
            dim = [geometry.Nx, geometry.Ny, geometry.Nz, self.Nlevel]

        self.n     = read_farray(dim, up, "double")
        self.nstar = read_farray(dim, up, "double")

        up.done()

    def readrates(self, geometry, filename):
        
        f  = open(filename, 'rb')
        up = xdrlib.Unpacker(f.read())
        f.close()

        Nrad = self.Nline + self.Ncont
        
        if geometry.type == "ONE_D_PLANE":   
            dim = [geometry.Ndep]
        elif geometry.type == 'SPHERICAL_SYMMETRIC':
            dim = [geometry.Nradius]
        elif geometry.type == 'TWO_D_PLANE':
            dim = [geometry.Nx, geometry.Nz]
        elif geometry.type == 'THREE_D_PLANE':
            dim = [geometry.Nx, geometry.Ny, geometry.Nz]

        for kr in range(Nrad):
            self.transition[kr].Rij = read_farray(dim, up, "double")
            self.transition[kr].Rji = read_farray(dim, up, "double")
        up.done()
 
    def readdamping(self, geometry, filename):
        
        f  = open(filename, 'rb')
        up = xdrlib.Unpacker(f.read())
        f.close()
         
        if geometry.type == "ONE_D_PLANE":   
            dim = [geometry.Ndep]
        elif geometry.type == 'SPHERICAL_SYMMETRIC':
            dim = [geometry.Nradius]
        elif geometry.type == 'TWO_D_PLANE':
            dim = [geometry.Nx, geometry.Nz]
        elif geometry.type == 'THREE_D_PLANE':
            dim = [geometry.Nx, geometry.Ny, geometry.Nz]

        self.vbroad = read_farray(dim, up, "double")
        
        for kr in range(self.Nline):
            self.transition[kr].adamp =\
                read_farray(dim, up, "double")
        up.done()
        
    def readcollisions(self, geometry, filename):
        
        f  = open(filename, 'rb')
        up = xdrlib.Unpacker(f.read())
        f.close()

        if geometry.type == "ONE_D_PLANE":   
            dim = [geometry.Ndep, self.Nlevel, self.Nlevel]
        elif geometry.type == 'SPHERICAL_SYMMETRIC':
            dim = [geometry.Nradius, self.Nlevel, self.Nlevel]
        elif geometry.type == 'TWO_D_PLANE':
            dim = [geometry.Nx, geometry.Nz, self.Nlevel, self.Nlevel]
        elif geometry.type == 'THREE_D_PLANE':
            dim = [geometry.Nx, geometry.Ny, geometry.Nz,\
                   self.Nlevel, self.Nlevel]

        self.Cij = read_farray(dim, up, "double")
        
        up.done()
 
    class transition:

        def __init__(self, up):
            self.read(up)

        def read(self, up):

            shapes = {0:"GAUSS", 1: "VOIGT", 2: "PRD",\
                      3: "HYDROGENIC", 4: "EXPLICIT"}
        
            types  = {0: "ATOMIC_LINE", 1: "ATOMIC_CONTINUUM"}
        
            self.type     = types[up.unpack_int()]
            self.i        = up.unpack_int()
            self.j        = up.unpack_int()
            self.Nwave    = up.unpack_int()
            self.Nblue    = up.unpack_int()
            self.lambda0  = up.unpack_double()            
            self.shape    = shapes[up.unpack_int()]
            self.strength = up.unpack_double()

            
    class fixed:
        def __init__(self, up):
            self.read(up)

        def read(self, up):
            self.type     = up.unpack_int()
            self.option   = up.unpack_int()
            self.i        = up.unpack_int()
            self.j        = up.unpack_int()
            self.lambda0  = up.unpack_double()
            self.strength = up.unpack_double()
            self.Trad     = up.unpack_double()

