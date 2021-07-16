import xdrlib
import numpy as np
import os

import rhanalyze.rhgeometry
from rhanalyze.rhtools import read_farray, read_string

class backgrflags:

    def __init__(self, hasline, ispolarized):
        self.hasline = hasline
        self.ispolarized = ispolarized

class opacity:
    
    def __init__(self, inputs, geometry, atmos, spectrum, rhdir):
        
        self.rhdir    = rhdir
        self.inputs   = inputs
        self.geometry = geometry
        self.atmos    = atmos
        self.spectrum = spectrum

        self.read_BRS()
        self.read_ASRS()

    def read_BRS(self):

        brsfile = '{0}/brs.out'.format(self.rhdir)
        f  = open(brsfile, 'rb')
        up = xdrlib.Unpacker(f.read())
        f.close()

        atmosID = read_string(up)
        Nspace  = up.unpack_int()
        
        Nspect  = up.unpack_int()

        hasline     = read_farray(Nspect, up, "int")
        ispolarized = read_farray(Nspect, up, "int")
        self.backgrflags = backgrflags(hasline, ispolarized)
        
        if self.atmos.moving or self.atmos.stokes:
            dim = 2*Nspect * self.geometry.Nrays
        else:
            dim = Nspect

        self.bg_recno = read_farray(dim, up, "int")
        up.done()
        
    def read_ASRS(self):

        asrsfile = '{0}/asrs.out'.format(self.rhdir)
        if not os.path.isfile(asrsfile):
            self.as_recno = None
            return
        
        f  = open(asrsfile, 'rb')
        up = xdrlib.Unpacker(f.read())
        f.close()

        if self.atmos.moving or self.atmos.stokes or\
           self.inputs.PRD_angle_dep:
            dim = self.spectrum.Nspect * self.geometry.Nrays
        else:
            dim = self.spectrum.Nspect

        self.as_recno = read_farray(dim, up, "int")
        up.done()
                    
    def read(self, waveno, rayno):

        SIZE_OF_DOUBLE = 8

        if waveno >= self.spectrum.Nspect:
            print("waveno {0} >= "
                  "spectrum.Nspect = {1}".format(waveno, self.spectrum.Nspect))
            waveno = 0
            
        if rayno >= self.geometry.Nrays:
            print("rayno {0} >= "
                  "geometry.Nrays = {1}".format(rayno, self.geometry.Nrays))
            rayno = 0
            
        self.waveno = waveno
        self.rayno  = rayno

        if self.geometry.type == "ONE_D_PLANE":           
            reclen_as = 2 * self.geometry.Ndep * SIZE_OF_DOUBLE
            reclen_bg = self.geometry.Ndep * SIZE_OF_DOUBLE
            
            dim1 = self.geometry.Ndep
            
        elif self.geometry.type == 'SPHERICAL_SYMMETRIC':
            reclen_as = 2 * self.geometry.Nradius * SIZE_OF_DOUBLE
            reclen_bg = self.geometry.Nradius * SIZE_OF_DOUBLE

            dim1 = self.geometry.Nradius
            
        elif self.geometry.type == 'TWO_D_PLANE':
            reclen_as = 2 * (self.geometry.Nx*self.geometry.Nz) * SIZE_OF_DOUBLE
            reclen_bg = (self.geometry.Nx*self.geometry.Nz) * SIZE_OF_DOUBLE

            dim1 = self.geometry.Nx * self.geometry.Nz
            
        elif self.geometry.type == 'THREE_D_PLANE':
            reclen_as = 2 * (self.geometry.Nx*self.geometry.Ny*\
                             self.geometry.Nz) * SIZE_OF_DOUBLE
            reclen_bg = (self.geometry.Nx*self.geometry.Ny*\
                         self.geometry.Nz) * SIZE_OF_DOUBLE

            dim1 = self.geometry.Nx * self.geometry.Ny *self.geometry.Nz
            

        if self.atmos.moving or self.atmos.stokes or \
           self.inputs.PRD_angle_dep:
            index_as = waveno*self.geometry.Nrays + rayno
        else:
            index_as = waveno

        if self.atmos.moving or self.atmos.stokes:
            index_bg = 2 * (waveno * self.geometry.Nrays + rayno) + 1
        else:
            index_bg = waveno

            
        file_as = open('{0}/opacity.out'.format(self.rhdir), 'rb')

        offset_as = self.as_recno[index_as] * reclen_as
        chunk_as  = 2 * dim1 * SIZE_OF_DOUBLE
        file_as.seek(offset_as, 0)
        up_as = xdrlib.Unpacker(file_as.read(chunk_as))

        self.chi_as = read_farray(dim1, up_as, "double")
        self.eta_as = read_farray(dim1, up_as, "double")
        file_as.close()
        up_as.done()

        
        file_bg = open('{0}/background.dat'.format(self.rhdir), 'rb')
        offset_bg = self.bg_recno[index_bg] * reclen_bg
        file_bg.seek(offset_bg, 0)
        
        if self.atmos.stokes and self.backgrflags.ispolarized[waveno]:
            chunk_bg = 9 * dim1 * SIZE_OF_DOUBLE
            if self.inputs.magneto_optical:
                chunk_bg += 3 * dim1 * SIZE_OF_DOUBLE
        else:
            chunk_bg = 3 * dim1 * SIZE_OF_DOUBLE

        buffer_bg = file_bg.read(chunk_bg)   
        self.chi_c = np.frombuffer(buffer_bg, dtype='float', count=dim1)

        if self.atmos.stokes and self.backgrflags.ispolarized[waveno]:
            offset = 4 * dim1 * SIZE_OF_DOUBLE
            if self.inputs.magneto_optical:
                offset += 3 * dim1 * SIZE_OF_DOUBLE
        else:
            offset = dim1 * SIZE_OF_DOUBLE
            
        self.eta_c = np.frombuffer(buffer_bg, dtype='float', count=dim1, \
                                   offset=offset)
        
        if self.atmos.stokes and self.backgrflags.ispolarized[waveno]:
            offset += 3 * dim1 * SIZE_OF_DOUBLE
        else:
            offset += dim1 * SIZE_OF_DOUBLE

        self.scatt = np.frombuffer(buffer_bg, dtype='float', count=dim1, \
                                   offset=offset)


        if self.geometry.type == 'TWO_D_PLANE':
            np.reshape(self.chi_as, [self.geometry.Nx, self.geometry.Nz])
            np.reshape(self.eta_as, [self.geometry.Nx, self.geometry.Nz])

            np.reshape(self.chi_c, [self.geometry.Nx, self.geometry.Nz])
            np.reshape(self.eta_c, [self.geometry.Nx, self.geometry.Nz])
            np.reshape(self.scatt, [self.geometry.Nx, self.geometry.Nz])
            
        elif self.geometry.type == 'THREE_D_PLANE':
            np.reshape(self.chi_as, [self.geometry.Nx, self.geometry.Ny,\
                                     self.geometry.Nz])
            np.reshape(self.eta_as, [self.geometry.Nx, self.geometry.Ny,\
                                     self.geometry.Nz])

            np.reshape(self.chi_c, [self.geometry.Nx, self.geometry.Ny,\
                                    self.geometry.Nz])
            np.reshape(self.eta_c, [self.geometry.Nx, self.geometry.Ny,\
                                    self.geometry.Nz])
            np.reshape(self.scatt, [self.geometry.Nx, self.geometry.Ny,\
                                    self.geometry.Nz])
            
        
        self.read_J()

    def read_J(self):

        SIZE_OF_DOUBLE = 8

        if self.geometry.type == "ONE_D_PLANE":           
            dim1 = self.geometry.Ndep
        elif self.geometry.type == 'SPHERICAL_SYMMETRIC':
            dim1 = self.geometry.Nradius
        elif self.geometry.type == 'TWO_D_PLANE':
            dim1 = self.geometry.Nx * self.geometry.Nz
        elif self.geometry.type == 'THREE_D_PLANE':
            dim1 = self.geometry.Nx * self.geometry.Ny *self.geometry.Nz

        Jfile  = '{0}/J.dat'.format(self.rhdir)
        offset = self.waveno * dim1 * SIZE_OF_DOUBLE
        self.J = np.fromfile(Jfile, dtype='float', count=dim1, offset=offset)

        if self.geometry.type == 'TWO_D_PLANE':
            np.reshape(self.J, [self.geometry.Nx, self.geometry.Nz])
        elif self.geometry.type == 'THREE_D_PLANE':
            np.reshape(self.J, [self.geometry.Nx, self.geometry.Ny,\
                                self.geometry.Nz])
            
    def Source(self):

        chi_tot = self.chi_as + self.chi_c
        self.S  = (self.eta_as + self.eta_c + self.J * self.scatt) / chi_tot

        
    def Planck(self, HZ=True):
        
        CLIGHT     = 2.99792458E+08
        HPLANCK    = 6.626176E-34
        KBOLTZMANN = 1.380662E-23

        CM_TO_M      = 1.0E-02
        NM_TO_M      = 1.0E-09
        ERG_TO_JOULE = 1.0E-07

        lambda_m   = NM_TO_M * self.spectrum.waves[self.waveno]

        hc_kl      = (HPLANCK * CLIGHT) / (KBOLTZMANN * lambda_m)
        twohnu3_c2 = (2.0*HPLANCK*CLIGHT) / lambda_m**3

        if HZ:
            
            ## --- In [J m^-2 s^-1 Hz^-1 sr^-1] --     --------------- ;

            self.Bp = twohnu3_c2 / (np.exp(hc_kl/self.atmos.T) - 1.0)
        else:

            ## --- In [erg cm^-2 s^-1 nm^-1 sr^-1] --  --------------- ;

            C = NM_TO_M * CM_TO_M^2 / ERG_TO_JOULE

            self.Bp =  C * twohnu3_c2 * (CLIGHT/lambda_m**2) / \
                (np.exp(hc_kl/self.atmos.T) - 1.0)

    def get_tau(self, center=False):

        if center:
            xmu = 1.0
        else:
            xmu = self.geometry.xmu[rayno]

        if self.geometry.type == "ONE_D_PLANE":
            N    = self.geometry.Ndep
            tau  = np.zeros(N, dtype='float')
            path = self.geometry.height / xmu

        chi = self.chi_c + self.chi_as
            
        tau[0] = 0.0
        for k in range(1, N, 1):
            dtau = 0.5*(chi[k-1] + chi[k]) * (path[k-1] - path[k])
            tau[k] = tau[k-1] + dtau
            
        self.tau = tau
