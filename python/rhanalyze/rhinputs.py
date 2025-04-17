import sys

if sys.version_info[1] < 13:
    import xdrlib
else:
    import mda_xdrlib.xdrlib as xdrlib

    
class inputs:
    
    def __init__(self, filename='input.out'):
        self.filename = filename
        self.read()
         
    def read(self):
        
        stokes_mode = ["NO_STOKES", "FIELD_FREE", "FULL_STOKES"]
        
        f  = open(self.filename, 'rb')
        up = xdrlib.Unpacker(f.read())
        f.close()
        
        self.magneto_optical = up.unpack_int()
        self.PRD_angle_dep   = up.unpack_int()
        self.XRD             = up.unpack_int()
        self.start_solution  = up.unpack_int()
        
        sm = up.unpack_int()
        self.stokes_mode     = stokes_mode[sm]
        
        self.metallicity     = float(up.unpack_double())
        self.backgr_pol      = up.unpack_int()
        self.big_endian      = up.unpack_int()
        
        up.done()
    
