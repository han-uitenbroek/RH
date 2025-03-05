import sys

if sys.version_info[1] < 13:
    import xdrlib
else:
    import mda_xdrlib

import numpy as np

from rhanalyze.rhtools import read_farray


class geometry:
    
    def __init__(self, filename='geometry.out'):
        self.filename = filename
        self.read()
         
    def read(self):

        geometry_types = ['ONE_D_PLANE', 'TWO_D_PLANE',
                          'SPHERICAL_SYMMETRIC', 'THREE_D_PLANE']

        f  = open(self.filename, 'rb')
        up = xdrlib.Unpacker(f.read())
        f.close()

        self.type  = geometry_types[up.unpack_int()]
        self.Nrays = up.unpack_int()

        if self.type == 'ONE_D_PLANE':

            self.Ndep = up.unpack_int()

            self.xmu    = read_farray(self.Nrays, up, "double")
            self.wmu    = read_farray(self.Nrays, up, "double")
            
            self.height = read_farray(self.Ndep, up, "double")
            self.cmass  = read_farray(self.Ndep, up, "double")
            self.tau500 = read_farray(self.Ndep, up, "double")
            self.vz     = read_farray(self.Ndep, up, "double")
            
        elif self.type == 'SPHERICAL_SYMMETRIC':

            self.Nradius = up.unpack_int()
            self.Ncore   = up.unpack_int()
            self.radius  = up.unpack_double()

            self.xmu    = read_farray(self.Nrays, up, "double")
            self.wmu    = read_farray(self.Nrays, up, "double")
            
            self.r      = read_farray(self.Nradius, up, "double")
            self.cmass  = read_farray(self.Nradius, up, "double")
            self.tau500 = read_farray(self.Nradius, up, "double")
            self.vr     = read_farray(self.Nradius, up, "double")
            
        elif self.type == 'TWO_D_PLANE':

            self.Nx  = up.unpack_int()
            self.Nz  = up.unpack_int()

            self.AngleSet = up.unpack_int()
            self.xmu = read_farray(self.Nrays, up, "double")
            self.ymu = read_farray(self.Nrays, up, "double")
            self.wmu = read_farray(self.Nrays, up, "double")

            self.x = read_farray(self.Nx, up, "double")
            self.z = read_farray(self.Nz, up, "double")

            dim = [self.Nx, self.Nz]
            self.vx = read_farray(dim, up, "double")
            self.vz = read_farray(dim, up, "double")
            
        elif self.type == 'THREE_D_PLANE':
            
            self.Nx  = up.unpack_int()
            self.Ny  = up.unpack_int()
            self.Nz  = up.unpack_int()

            self.AngleSet = up.unpack_int()
            self.xmu = read_farray(self.Nrays, up, "double")
            self.ymu = read_farray(self.Nrays, up, "double")
            self.wmu = read_farray(self.Nrays, up, "double")

            self.dx = up.unpack_double()
            self.dy = up.unpack_double()
            self.z  = read_farray(self.Nz, up, "double")

            dim = [self.Nx, self.Ny, self.Nz]
            self.vx = read_farray(dim, up, "double")
            self.vy = read_farray(dim, up, "double")
            self.vz = read_farray(dim, up, "double")
        
        up.done()
    
