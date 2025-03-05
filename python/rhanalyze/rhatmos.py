import sys

if sys.version_info[1] < 13:
    import xdrlib
else:
    import mda_xdrlib

import numpy as np

import rhanalyze.rhgeometry
from rhanalyze.rhtools import read_farray, read_string, write_farray


class element:

    def __init__(self, up):
        self.read_element(up)
        
    def read_element(self, up):
        self.ID     = read_string(up)
        self.weight = float(up.unpack_double())
        self.abund  = float(up.unpack_double())
    
class atmos:
    
    def __init__(self, geometry, filename='atmos.out'):
        self.filename = filename
        self.read(geometry)
         
    def read(self, geometry):

        f  = open(self.filename, 'rb')
        up = xdrlib.Unpacker(f.read())
        f.close()

        self.NHydr  = up.unpack_int()
        self.Nelem  = up.unpack_int()
        self.moving = up.unpack_int()

        if geometry.type == "ONE_D_PLANE":           
            dim1 = [geometry.Ndep]         
            dim2 = [geometry.Ndep, self.NHydr]
            
        elif geometry.type == 'SPHERICAL_SYMMETRIC':
            dim1 = [geometry.Nradius]         
            dim2 = [geometry.Nradius, self.NHydr]
            
        elif geometry.type == 'TWO_D_PLANE':
            dim1 = [geometry.Nx, geometry.Nz]
            dim2 = [geometry.Nx, geometry.Nz, self.NHydr]

        elif geometry.type == 'THREE_D_PLANE':
            dim1 = [geometry.Nx, geometry.Ny, geometry.Nz]
            dim2 = [geometry.Nx, geometry.Ny, geometry.Nz, self.NHydr]


        self.T      = read_farray(dim1, up, "double")
        self.n_elec = read_farray(dim1, up, "double")
        self.vturb  = read_farray(dim1, up, "double")

        self.nH = read_farray(dim2, up, "double")       
        self.ID = read_string(up)

        
        self.elements = {}
        for n in range(self.Nelem):
            self.elements[n] = element(up)

        if geometry.type != 'SPHERICAL_SYMMETRIC':
            try:
                stokes = up.unpack_int()
            except EOFError or IOError:
                self.stokes = False
                return
            else:
                self.stokes = True

                self.B       = read_farray(dim1, up, "double")
                self.gamma_B = read_farray(dim1, up, "double")
                self.chi_B   = read_farray(dim1, up, "double")
                
        up.done()

class input_atmos:
    def __init__(self, geometrytype, atmosfile, Bfile=None, New=False):
        
        self.type      = geometrytype
        self.atmosfile = atmosfile
        self.Bfile     = Bfile
        
        if New == False:
            self.read()
    
    def read(self):
       
        if self.type == "ONE_D_PLANE" or self.type == "SPHERICAL_SYMMETRIC":

            CM_TO_M = 1.0E-2
            G_TO_KG = 1.0E-3

            data = []
            with open(self.atmosfile, 'r') as file:
                for line in file:
                    if line.startswith('*'):
                        continue
                    data.append(line.strip())
            
            self.ID    = data[0]
            scale      = data[1][0]
            
            if self.type == "ONE_D_PLANE":
                self.grav = float(data[2])
                self.Ndep = int(data[3])
                Nd = self.Ndep
            else:
                self.grav, self.radius = [float(x) for x in data[2].split()]
                self.Nradius, self.Ncore, self.Ninter = \
                    [int(x) for x in data[3].split()]
                Nd = self.Nradius

            self.grav  = np.power(10.0, self.grav) * CM_TO_M
            self.NHydr = 6

            hscale      = np.array(range(Nd), dtype="float")
            self.T      = np.array(range(Nd), dtype="float")
            self.n_elec = np.array(range(Nd), dtype="float")
            self.v      = np.array(range(Nd), dtype="float")
            self.vturb  = np.array(range(Nd), dtype="float")

            for n in range(Nd):
                hscale[n], self.T[n],\
                    self.n_elec[n], self.v[n], self.vturb[n] =\
                        [float(x) for x in data[n+4].split()]
            
            if scale == 'M':
                self.scale  = 'MASS_SCALE'
                self.cmass  = np.power(10.0, hscale)
                self.cmass *= G_TO_KG / CM_TO_M**2

            elif scale == 'T':
                self.scale  = 'TAU500_SCALE'
                self.tau500 = np.power(10.0, hscale)
                
            elif scale == 'H':
                self.scale  = 'GEOMETRIC_SCALE'
                self.height = hscale

            if len(data) > (4 + Nd):
                self.HLTE = False
                self.nH = np.array(range(Nd * self.NHydr),\
                                   dtype="float").reshape([Nd,\
                                                           self.NHydr],\
                                                          order='F')
                for n in range(Nd):
                    self.nH[n,:] =\
                        [float(x) for x in data[n+4+Nd].split()]
            else:
                self.HLTE = True
                                  
            self.nH     /= CM_TO_M**3
            self.n_elec /= CM_TO_M**3

            dim1 = [Nd]

        elif self.type == "TWO_D_PLANE" or self.type == "THREE_D_PLANE":
            
            f = open(self.atmosfile, 'rb')
            up = xdrlib.Unpacker(f.read())
            f.close()
             
            if self.type == "TWO_D_PLANE":
                self.Nx = up.unpack_int()
                self.Nz = up.unpack_int()
                self.NHydr = up.unpack_int()

                self.boundary = read_farray([3], up, "int")

                self.dx = read_farray([self.Nx], up, "double")
                self.z  = read_farray([self.Nz], up, "double")

                dim1 = [self.Nx, self.Nz]
                dim2 = [self.Nx, self.Nz, self.NHydr]
                
            elif self.type == "THREE_D_PLANE":
                self.Nx = up.unpack_int()
                self.Ny = up.unpack_int()
                self.Nz = up.unpack_int()
                self.NHydr = up.unpack_int()

                self.boundary = read_farray([2], up, "int")

                self.dx = up.unpack_double()
                self.dy = up.unpack_double()
                self.z  = read_farray([self.Nz], up, "double")

                dim1 = [self.Nx, self.Ny, self.Nz]
                dim2 = [self.Nx, self.Ny, self.Nz, self.NHydr]
                
            self.T      = read_farray(dim1, up, "double")
            self.n_elec = read_farray(dim1, up, "double")
            self.vturb  = read_farray(dim1, up, "double")
            self.vx     = read_farray(dim1, up, "double")
            
            if self.type == "THREE_D_PLANE":
                self.vy     = read_farray(dim1, up, "double")
                
            self.vz     = read_farray(dim1, up, "double")

            self.nH     = read_farray(dim2, up, "double")
                
            up.done()
            
        else:
            print("Not a valid input atmosphere type: {0}".format(self.type))
            return
            
        if self.Bfile != None:
            
            f = open(self.Bfile, 'rb')
            up = xdrlib.Unpacker(f.read())
            f.close()

            self.B     = read_farray(dim1, up, "double")
            self.gamma = read_farray(dim1, up, "double")
            self.chi   = read_farray(dim1, up, "double")
            
            up.done()

    def write(self, outfile=None, Bfile=None):

        if self.type == "ONE_D_PLANE" or self.type == "SPHERICAL_SYMMETRIC":

            CM_TO_M = 1.0E-2
            G_TO_KG = 1.0E-3

            nH     = self.nH.copy() * CM_TO_M**3
            n_elec = self.n_elec.copy() * CM_TO_M**3

            data = []

            data.append("* Model atmosphere written by " \
                        "rhatmos.input_atmos.write()\n")
            data.append("*\n")
            data.append("  {0}\n".format(self.ID))

            if self.scale == "MASS_SCALE":
                hscale = np.log10(self.cmass / (G_TO_KG / CM_TO_M**2))
                data.append("  Mass scale\n")
            elif self.scale == "TAU500_SCALE":
                hscale = np.log10(self.tau500)
                data.append("  Tau500 scale\n")
            elif self.scale == "GEOMETRIC_SCALE":
                hscale = self.height
                data.append("  Height scale\n")

            data.append('*\n')

            grav = np.log10(self.grav / CM_TO_M)
    
            if self.type == "ONE_D_PLANE":
                data.append("* lg g [cm s^-2]\n")
                data.append('     {:5.2f}\n'.format(grav))
                data.append("* Ndep\n")   
                data.append('   {:4d}\n'.format(self.Ndep))
                
                Nd = self.Ndep
            else:
                data.append("* lg g [cm s^-2]      Radius [km]\n")
                data.append('     {:5.2f}          '\
                            '{:7.2E}\n'.format(grav, self.radius))
                data.append("* Nradius   Ncore   Ninter\n")
                fmt = 3 * '    {:4d}' + "\n"
                data.append(fmt.format(self.Nradius, self.Ncore, self.Ninter))
                            
                Nd = self.Nradius

            data.append("*\n")
            
            if self.scale == "MASS_SCALE":
                data.append("*  lg column Mass   Temperature    "\
                            "Ne             V              Vturb\n")
            elif self.scale == "TAU500_SCALE":
                data.append("*lg tau500          Temperature    "\
                            "Ne             V              Vturb\n")
            elif self.scale == "GEOMETRIC_SCALE":
                data.append("*height [km]        Temperature    "\
                            "Ne             V              Vturb\n")

            fmt = '  {: 12.8E}' + 4 * '  {: 10.6E}' + "\n"
            for k in range(Nd):
                data.append(fmt.format(hscale[k], self.T[k], n_elec[k],\
                                       self.v[k], self.vturb[k]))

            data.append("*\n")
            
            if not self.HLTE:
                data.append("* NLTE Hydrogen populations\n")
                data.append("*  nh[1]        nh[2]        nh[3]        "\
                            "nh[4]        nh[5]        np\n")

                fmt = self.NHydr * '   {:8.4E}' + "\n"
                for k in range(Nd):
                    data.append(fmt.format(*nH[k, :]))
                    
            if outfile == None:
                f = open(self.atmosfile, 'w')
            else:
                f = open(outfile, 'w')
                
            for line in data:
                f.write(line)
            f.close()

             
        elif self.type == "TWO_D_PLANE" or self.type == "THREE_D_PLANE":

            pck = xdrlib.Packer()

            if self.type == "TWO_D_PLANE":
                write_farray(np.array([self.Nx, self.Nz, self.NHydr]),\
                             pck, "int")
                write_farray(np.array(self.boundary), pck, "int")
                write_farray(self.dx, pck, "double")
               
            elif self.type == "THREE_D_PLANE":
                write_farray(np.array([self.Nx, self.Ny,\
                                       self.Nz, self.NHydr]),\
                             pck, "int")
                write_farray(np.array(self.boundary), pck, "int")
                pck.pack_double(self.dx)
                pck.pack_double(self.dy)

            write_farray(self.z, pck, "double")
            write_farray(self.T, pck, "double")
            write_farray(self.n_elec, pck, "double")
            write_farray(self.vturb, pck, "double")
            write_farray(self.vx, pck, "double")

            if self.type == "THREE_D_PLANE":
                write_farray(self.vy, pck, "double")
            
            write_farray(self.vz, pck, "double")
            write_farray(self.nH, pck, "double")
 
            if outfile == None:
                f = open(self.atmosfile, 'wb')
            else:
                f = open(outfile, 'wb')
                
            f.write(pck.get_buffer())
            f.close()
            pck.reset()
            
        else:
            print("Not a valid input atmosphere type: {0}".format(self.type))
            return

        if self.Bfile != None:

            pck = xdrlib.Packer()

            write_farray(self.B, pck, "double")
            write_farray(self.gamma, pck, "double")
            write_farray(self.chi, pck, "double")

            f = open(self.Bfile, 'wb')
            f.write(pck.get_buffer())
            f.close()
            pck.reset()
