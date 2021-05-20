import xdrlib
import numpy as np


def read_farray(dim, up, dtype="double"):
    
    if dtype == "float":
        func = up.unpack_float
        dt   = "float"
    elif dtype == "double":
        func = up.unpack_double
        dt   = "float"
    elif dtype == "int":
        func = up.unpack_int
        dt   = "int"

    N = np.prod(dim)
    farray = np.array(up.unpack_farray(N, func), dtype=dt)

    return np.reshape(farray, dim, order='F')


def write_farray(fa, pck, dtype="double"):

    if dtype == "float":
        func = pck.pack_float
    elif dtype == "double":
        func = pck.pack_double
    elif dtype == "int":
        func = pck.pack_int

    N      = np.prod(fa.shape)
    farray = np.reshape(fa, N, order='F')
    
    pck.pack_farray(N, farray, func)
    
    
def read_string(up, size=0):
    """Compensate for IDL string quirks"""

    up.unpack_int()

    if size == 0:
        return str(up.unpack_string().strip(), 'utf-8')
    else:
        up.unpack_int()
        return str(up.unpack_fstring(size), 'utf-8')
