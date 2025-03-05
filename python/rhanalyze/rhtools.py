import sys

if sys.version_info[1] < 13:
    import xdrlib
else:
    import mda_xdrlib

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
        return str(up.unpack_string().strip(), encoding='utf-8')
    else:
        up.unpack_int()
        return str(up.unpack_fstring(size), encoding='utf-8')


def monotonic_increasing(x):
    dx = np.diff(x)
    return np.all(dx > 0)

def monotonic_decreasing(x):
    dx = np.diff(x)
    return np.all(dx < 0)

def monotonic(x):
    dx = np.diff(x)
    return np.all(dx < 0) or np.all(dx > 0)

def table_invert(table, values, mode=None):
    
    if np.ndim(table) == 0:
        print("Table cannot be a scalar!")
        return table
        
    t_array = np.array(table)

    tsize = np.size(t_array)
    if  tsize <= 1:
        print("Table cannot have one element!")
        return t_array
    
    if monotonic(t_array) == False:
        print("Table has to be monotonic!")
        return None

    increasing = True
    if t_array[0] > t_array[1]:
        increasing = False
        
    if np.ndim(values) == 0:
        v_array = np.array([values])
    else:  
        v_array = np.array(values)
    
    lookup = []
    for v in v_array:
        if increasing:
            
            if v <= t_array[0]:
                lookup.append({"index": 0, "value": t_array[0], "eff_index": 0})
            elif v >= t_array[-1]:
                lookup.append({"index": len(t_array)-1, "value": t_array[-1], "eff_index": len(t_array)-1})
            else:
                index = np.argmin(np.abs(t_array - v))

                if v > t_array[index]:
                    eff_index = index + (v - t_array[index]) / (t_array[index+1] - t_array[index])
                else:
                    eff_index = index - 1 + \
                    (v - t_array[index-1]) / (t_array[index] - t_array[index-1])
            
                lookup.append({"index": index, "value": v, "eff_index": eff_index})

        else:

            if v >= t_array[0]:
                lookup.append({"index": 0, "value": t_array[0], "eff_index": 0})
            elif v <= t_array[-1]:
                lookup.append({"index": len(t_array)-1, "value": t_array[-1], "eff_index": len(t_array)-1})
            else:
                index = np.argmin(np.abs(t_array - v))

                if v < t_array[index]:
                    eff_index = index + (v - t_array[index]) / (t_array[index+1] - t_array[index])
                else:
                    eff_index = index - 1 + \
                    (t_array[index-1] - v) / (t_array[index-1] - t_array[index])
            
                lookup.append({"index": index, "value": v, "eff_index": eff_index})
    
    if mode == "effective":
        return np.array([dict["eff_index"] for dict in lookup])
    elif mode == "value":
        return np.array([dict["value"] for dict in lookup])
    elif mode == "index":
        return np.array([dict["index"] for dict in lookup])
    else:
        return lookup
