"""
Main class definition for rhanalyze
"""
import os
import fnmatch

import rhanalyze.rhinputs
import rhanalyze.rhgeometry
import rhanalyze.rhatmos
import rhanalyze.rhspectrum
import rhanalyze.rhatom
import rhanalyze.rhopacity


class rhout:

    def __init__(self, rhdir='.'):

        self.rhdir = rhdir
        
        self.inputs = rhinputs.inputs(os.path.join(self.rhdir, "input.out"))

        self.geometry = rhgeometry.geometry(os.path.join(self.rhdir, "geometry.out"))
        self.atmos = rhatmos.atmos(self.geometry, \
                                   filename=os.path.join(self.rhdir, "atmos.out"),
                                   molfile=os.path.join(self.rhdir, "molecules.out"))
        
        self.spectrum = rhspectrum.spectrum(self.inputs,\
                                            self.geometry,\
                                            self.atmos,\
                                            filename=os.path.join(self.rhdir, "spectrum.out"), \
                                            fluxfile=os.path.join(self.rhdir, "flux.out"))

        self.rays = []
        for file in sorted(os.listdir(rhdir)):
            if fnmatch.fnmatch(file, 'spectrum_?.??*'):
                self.rays.append(rhspectrum.rays(self.inputs,\
                                                 self.geometry,\
                                                 self.spectrum,\
                                                 filename=os.path.join(self.rhdir, file)))
        self.Nray = len(self.rays)

        self.atoms = []
        for file in sorted(os.listdir(rhdir)):
            if fnmatch.fnmatch(file, 'atom.*.out'):
                self.atoms.append(rhatom.atoms(self.geometry, \
                                               path=os.path.join(self.rhdir, file)))
        self.Natom = len(self.atoms)

        self.opacity = rhopacity.opacity(self.inputs, self.geometry,\
                                         self.atmos, self.spectrum, \
                                         self.rhdir)
