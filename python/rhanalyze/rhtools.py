import sys
import os
import glob
import numpy as np


if sys.version_info[1] < 13:
    import xdrlib
else:
    import mda_xdrlib.xdrlib as xdrlib


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


### Class and methods to read, write and manipulate RH-style atomic models

class RHatom:
    
    def __init__(self, file=None):
        
        self.atom_ID = ""
        self.levels  = []
        
        if (file != None):
            self.read(file)

        
    def read(self, atom_file):
        
        ## Read the file
        self.atom_file = atom_file
        file = open(self.atom_file, "r")

        ## Strip all empty and comment lines
        lines = [line for line in file.readlines() if (line.strip() and not line.startswith("#"))]
        file.close()

        lineNo = 0
        self.atom_ID = lines[lineNo].strip()
        self.Nlevel, self.Nbb, self.Nbf, self.Nfix = np.fromstring(lines[1], dtype=int, sep=' ')
        lineNo += 2

        for i in range(self.Nlevel):
            self.levels.append(self.level())
            self.levels[i].parse(lines[lineNo])
            lineNo += 1

        self.levels.sort(key=lambda level: level.E)
        labels = [level.label for level in self.levels]
        
        self.bb_trans = []
        for l in range(self.Nbb):

            bb = self.boundbound()
            (i, j) = bb.parse(lines[lineNo])
            
            bb.label_i = labels[i]
            bb.label_j = labels[j]
            if self.levels[i].E > self.levels[j].E:
                (bb.label_i, bb.label_j) = (bb.label_j, bb.label_i)

            self.bb_trans.append(bb)
            lineNo += 1

        self.fix_trans = []
        for l in range(self.Nfix):

            fix = self.fixed()
            (i, j) = fix.parse(lines[lineNo])
            
            fix.label_i = labels[i]
            fix.label_j = labels[j]
            if self.levels[i].E > self.levels[j].E:
                (fix.label_i, fix.label_j) = (fix.label_j, fix.label_i)

            self.fix_trans.append(fix)
            lineNo += 1

        self.bf_trans = []
        for l in range(self.Nbf):

            bf = self.boundfree()
            (i, j, lines_read) = bf.parse(lines[lineNo:])
            
            bf.label_i = labels[i]
            bf.label_j = labels[j]
            if self.levels[i].E > self.levels[j].E:
                (bf.label_i, bf.label_j) = (bf.label_j, bf.label_i)

            self.bf_trans.append(bf)
            lineNo += lines_read

        end = False
        T   = None
        self.coll_trans = []
        
        while end == False:
            xsection = self.coll_cross()
            if not np.isscalar(T):
                xsection.T = T
            (i, j, end) = xsection.parse(lines[lineNo])

            if end == True:
                return
                
            if xsection.entry == "TEMP":
                T = xsection.T
            else:
                xsection.label_i = labels[i]
                xsection.label_j = labels[j]
                if self.levels[i].E > self.levels[j].E:
                    (xsection.label_i, xsection.label_j) = (xsection.label_j, xsection.label_i)

                xsection.T = T
                self.coll_trans.append(xsection)

            lineNo += 1
            
    def write(self, out_file, overwrite=False):

        if overwrite:
            file = open(out_file, "w")
        else:       
            try:
                file = open(out_file, "x")
            except FileExistsError:
                print("File already exists: ", out_file)
            return
            
        labels = [level.label for level in self.levels]
 
        lines = []
        lines.append("  {0:2s}\n".format(self.atom_ID))

        fmt = 4 * "    {:4d}" + '\n'
        lines.append(fmt.format(self.Nlevel, self.Nbb, self.Nbf, self.Nfix))

        fmt = " {0:>9.3f}   {1:>5.2f}    '{2:20s}'    {3:>2d}\n"
        for level in self.levels:
            lines.append(fmt.format(level.E, level.g, level.label, level.stage))

        fmt = 2*"{:3d}" + "   {:9.3E}  {:>6s}  {:>4d}  {:>6s}  {:>4.1f}  {:>6.1f}  {:>6s} " + \
              4*" {:5.2f}" + 2*"   {:7.2E}" + "\n"
        for bb in self.bb_trans:

            (i, j) = (labels.index(bb.label_i), labels.index(bb.label_j))
            lines.append(fmt.format(j, i, bb.log_gf, bb.line_shape, bb.Nlambda, bb.symmetr, bb.qcore, \
                         bb.qwing, bb.vdWaalsapprx, *(bb.vdWaals), bb.gamma_rad, bb.gamma_Stark))

        fmt = 2*"{:3d}" + "   {:9.3E}      {:>7.1f}       {:>20s}" + "\n"
        for fix in self.fix_trans:
            
            (i, j) = (labels.index(fix.label_i), labels.index(fix.label_j))
            lines.append(fmt.format(j, i, fix.strength, fix.Trad, fix.option))

        fmt1 = 2 * " {:>3d}" + "   {:10.4E}    {:>3d}     {:s}      {:6.1f}" + "\n"
        fmt2 = "   {:6.1f}    {:10.4E}" + "\n"
        for bf in self.bf_trans:

            (i, j) = (labels.index(bf.label_i), labels.index(bf.label_j))
            lines.append(fmt1.format(j, i, bf.alpha_0, bf.Nlambda, bf.bf_type, bf.lamb_min))

            if bf.bf_type == "EXPLICIT":
                for n in range(bf.Nlambda):
                    lines.append(fmt2.format(bf.waves[n], bf.alpha[n]))

        T = np.array([0.0])
        entry = ""
        for coll in self.coll_trans:

            (i, j) = (labels.index(coll.label_i), labels.index(coll.label_j))
            if (coll.T != T).all() or coll.entry != entry:
                fmt = " {:>7s}     {:3d}    " + len(coll.T) * "     {:9.1f}" + "\n"
                lines.append(fmt.format("TEMP", len(coll.T), *coll.T))
                T = coll.T

            fmt = " {:>7s}" + 2 * "  {:3d}" + len(coll.T) * "     {:9.3E}" + "\n"
            lines.append(fmt.format(coll.entry, j, i, *coll.xsection))
            entry = coll.entry

        lines.append("END\n")
        for line in lines:
            file.write(line)
        file.close()

        
    def add_level(self, label, E, g, stage):

        Elev = [level.E for level in self.levels]
        if E in Elev:
            print("Energy level ", E, " already exists")
            return
            
        labels = [level.label for level in self.levels]
        if label not in labels:
            level = self.level(E=E, label=label, g=g, stage=stage)

            self.levels.append(level)
            self.levels.sort(key=lambda level: level.E)
            self.Nlevel = len(self.levels)
        else:
            print("Label ", label, " already exists")
            return

    def delete_level(self, label):

        labels = [level.label for level in self.levels]
        
        if label not in labels:
            print("Label ", label, " does not exist")
        else:
            self.levels.pop(labels.index(label))
            self.Nlevel = len(self.levels)

        for bb in self.bb_trans:
            if bb.label_i == label or bb.label_j == label:
                self.bb_trans.remove(bb)
        self.Nbb = len(self.bb_trans)

        for bf in self.bf_trans:
            if bf.label_i == label or bf.label_j == label:
                self.bf_trans.remove(bf)
        self.Nbf = len(self.bf_trans)

        for coll in self.coll_trans:
            if coll.label_i == label or coll.label_j == label:
                self.coll_trans.remove(coll)

    def delete_bb(self, label_i, label_j):
        Ndel = 0
        for bb in self.bb_trans:
            if ((bb.label_i == label_i  and  bb.label_j == label_j) or \
                (bb.label_j == label_i  and  bb.label_i == label_j)):
                self.bb_trans.remove(bb)
                Ndel = 1
                break
        if Ndel == 0:
            print("No transition matched request: ", label_i, label_j)
        else:
            self.Nbb = len(self.bb_trans)

    def delete_bf(self, label_i, label_j):
        Ndel = 0
        for bf in self.bf_trans:
            if ((bf.label_i == label_i  and  bf.label_j == label_j) or \
                (bf.label_j == label_i  and  bf.label_i == label_j)):
                self.bf_trans.remove(bf)
                Ndel = 1
                break
        if Ndel == 0:
            print("No transition matched request: ", label_i, label_j)
        else:
            self.Nbb = len(self.bb_trans)

    def add_bb(self, label_i, label_j, log_gf, line_shape, Nlambda=15, symmetr="ASYMM", qcore=10, qwing=100, \
               vdwaalsapprx="UNSOLD", vdWaals=np.array([1.0, 0.0, 1.0, 0.0]), gamma_rad=1.0E8, gamma_Stark=1.0):
        labels = [level.label for level in self.levels]

        if label_i not in labels:
            print("add_bb: Label ", label_i, " not in level list")
            return
        if label_j not in labels:
            print("add_bb: Label ", label_j, " not in level list")
            return

        (i, j) = (labels.index(label_i), labels.index(label_j))
        if self.levels[i].E > self.levels[j].E:
            (label_i, label_j) = (label_j, label_i)
        pairs = [(bb.label_j, bb.label_i) for bb in self.bb_trans]
        if (label_j, label_i) in pairs:
            print("Bound-bound transition ", (label_j, label_i), " already exists")
            return
        
        bb = self.boundbound(label_i=label_i, label_j=label_j, log_gf=log_gf, line_shape=line_shape, \
                             Nlambda=Nlambda, \
                             symmetr=symmetr, qcore=qcore, qwing=qwing, \
                             vdWaalsapprx=vdWaalsapprx, vdWaals=vdWaals, \
                             gamma_rad=gamma_rad, gamma_Stark=gamma_Stark)
                             
        self.bb_trans.append(bb)
        self.Nbb = len(self.bb_trans)
        

    def add_bf(self, label_i, label_j, log_gf, alpha_0, bf_type, Nlambda, lambda_min):
        labels = [level.label for level in self.levels]

        if label_i not in labels:
            print("add_bf: Label ", label_i, " not in level list")
            return
        if label_j not in labels:
            print("add_bf: Label ", label_j, " not in level list")
            return

        (i, j) = (labels.index(label_i), labels.index(label_j))
        if self.levels[i].E > self.levels[j].E:
            (index_i, index_j) = (index_j, index_i)
        pairs = [(bf.label_j, bf.label_i) for bf in self.bf_trans]
        if (label_j, label_i) in pairs:
            print("Bound-free transition ", (label_j, label_i), " already exists")
            return

        bf = self.boundfree(label_i=label_i, label_j=label_j, alpha_0=alpha_0, bf_type=bf_type, \
                            Nlambda=Nlambda, lambda_min=lambda_min)
        self.bf_trans.append(bf)
        self.Nbf = len(self.bf_trans)


    class level:
        def __init__(self, E=None, g=None, label=None, stage=None):
            
            self.E = E
            self.g = g
            self.label = label
            self.stage = stage

        def parse(self, line):

            parts = line.split("'")
            self.E, self.g = np.fromstring(parts[0], dtype=float, sep=' ')
            self.label = parts[1]
            self.stage, levelNo = np.fromstring(parts[2], dtype=int, sep=' ')

    class boundbound:
        def __init__(self, label_i=None, label_j=None, log_gf=None, line_shape=None, \
                     Nlambda=15, symmetr="ASYMM", qcore=10, qwing=100, \
                     vdWaalsapprx="UNSOLD", vdWaals=np.array([1.0, 0.0, 1.0, 0.0]), \
                     gamma_rad=1.0E8, gamma_Stark=1.0):
            
            self.label_i = label_i
            self.label_j = label_j
            self.log_gf = log_gf
            self.line_shape = line_shape
            self.Nlambda = Nlambda
            self.symmetr = symmetr
            self.qcore = qcore
            self.qwing = qwing
            self.vdWaalsapprx = vdWaalsapprx
            self.vdWaals = vdWaals
            self.gamma_rad = gamma_rad
            self.gamma_Stark = gamma_Stark

        def parse(self, line):
        
            parts = line.split()
            
            i = np.int16(parts[0])
            j = np.int16(parts[1])
            self.log_gf = np.float64(parts[2])
            self.line_shape = parts[3]

            self.Nlambda = np.int16(parts[4])
            self.symmetr = parts[5]
            self.qcore   = np.float32(parts[6])
            self.qwing   = np.float32(parts[7])

            self.vdWaalsapprx = parts[8]
            self.vdWaals      = [np.float32(parts[9+n]) for n in range(4)]
            self.gamma_rad    = np.float32(parts[13])
            self.gamma_Stark  = np.float32(parts[14])

            if len(parts) > 15:
                self.gLande = np.float32(parts[15])

            return (i, j)

    class boundfree:

        def __init__(self, label_i=None, label_j=None, alpha_0=None, bf_type=None, \
                            Nlambda=None, lambda_min=None):

            self.label_i = label_i
            self.label_j = label_j
            self.alpha_0 = alpha_0

        def parse(self, lines):

            lines_read = 0
            parts = lines[lines_read].split()

            i = np.int16(parts[0])
            j = np.int16(parts[1])
            self.alpha_0 = np.float64(parts[2])
            self.Nlambda = np.int16(parts[3])
            self.bf_type = parts[4]
            self.lamb_min = np.float64(parts[5])
            lines_read += 1

            if (self.bf_type == 'EXPLICIT'):
                self.waves = np.zeros(self.Nlambda, dtype=np.float64)
                self.alpha = np.zeros(self.Nlambda, dtype=np.float64)
                for n in range(self.Nlambda):
                    parts = lines[lines_read].split()
                    self.waves[n] = np.float64(parts[0])
                    self.alpha[n] = np.float64(parts[1])
                    lines_read += 1
                        
            return(i, j, lines_read)

    class coll_cross:

        def __init__(self):
            
            self.label_i = None
            self.label_j = None
            self.T = None
            self.xsection = None

        def parse(self, line):

            end = False
            
            parts = line.split()
            self.entry = parts[0]

            (i, j) = (-1, -1)
            match self.entry:
                case "TEMP":
                    Nentry = np.int16(parts[1])
                    T = np.zeros(Nentry, dtype=np.float32)
                    for n in range(Nentry):
                        T[n] = np.float32(parts[2+n])
                    self.T = T
                    
                case "OMEGA" | "CI" | "CE":
                    i = np.int16(parts[1])
                    j = np.int16(parts[2])

                    self.xsection = np.zeros(len(self.T), dtype=np.float32)
                    for n in range(len(self.T)):
                        self.xsection[n] = np.float32(parts[3+n])

                case "END":
                    end = True
                case _:
                    print("Entry: " + self.entry + " not yet implemented")

            return (i, j, end)

    class fixed:

        def __init__(self, label_i=None, label_j=None, strength=None, Trad=None, option=None):

            self.label_i  = label_i
            self.label_j  = label_j
            self.strength = strength
            self.Trad     = Trad
            self.option   = option

        
        def parse(self, line):
            
            parts = line.split()

            i = np.int16(parts[0])
            j = np.int16(parts[1])
            self.strength = np.float64(parts[2])
            self.Trad     = np.float32(parts[3])
            self.option   = parts[4]

            return (i, j)


def vacuum_to_air(lambda_vac):

    #-# Source:
    #-#    [1] http://www.astro.uu.se/valdwiki/Air-to-vacuum%20conversion
    #-# Original reference:
    #-#    Donald Morton (2000, ApJ. Suppl., 130, 403)
    
    
    NM_TO_ANGSTROM = 10.0
    
    s2 = (1.0E4 / (NM_TO_ANGSTROM * lambda_vac))**2
    n  = 1.0 + 8.34254E-5 + 2.406147E-2 / (130.0 - s2) + 1.5998E-4 / (38.9 - s2)

    lambda_air = lambda_vac / n
    return lambda_air

def air_to_vacuum(lambda_air):

    #-# Source:
    #-#    [1] http://www.astro.uu.se/valdwiki/Air-to-vacuum%20conversion
    #-# Original reference:
    #-#    Donald Morton (2000, ApJ. Suppl., 130, 403)
    
    NM_TO_ANGSTROM = 10.0
    
    s2 = (1.0E4 / (NM_TO_ANGSTROM * lambda_air))**2
    n  = 1.0 + 8.336624212083E-5 + 2.408926869968E-2 / (130.1065924522 - s2) + \
        1.599740894897E-4 / (38.92568793293 - s2)
    
    lambda_vac = lambda_air * n
    return lambda_vac

    
class wavetable:

    def __init__(self, wavelengths):

        self.Nspect      = len(wavelengths)
        self.wavelengths = wavelengths

    def write(self, wavefile='wavetable.wave'):

        packer = xdrlib.Packer()

        packer.pack_int(self.Nspect)
        write_farray(air_to_vacuum(self.wavelengths), packer, dtype="double")

        f = open(wavefile, 'wb')
        f.write(packer.get_buffer())
        f.close()
