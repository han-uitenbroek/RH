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
        self.inputs = rhinputs.inputs('{0}/input.out'.format(rhdir))

        self.geometry = rhgeometry.geometry('{0}/geometry.out'.format(rhdir))
        self.atmos = rhatmos.atmos(self.geometry, '{0}/atmos.out'.format(rhdir))
        
        self.spectrum = rhspectrum.spectrum(self.inputs,\
                                            self.geometry,\
                                            self.atmos,\
                                            '{0}/spectrum.out'.format(rhdir))

        self.rays = []
        for file in sorted(os.listdir(rhdir)):
            if fnmatch.fnmatch(file, 'spectrum_?.??*'):
                self.rays.append(rhspectrum.rays(self.inputs,\
                                                 self.geometry,\
                                                 self.spectrum,\
                                                 filename=('{0}/'+file).format(rhdir)))
        self.Nray = len(self.rays)

        self.atoms = []
        for file in sorted(os.listdir(rhdir)):
            if fnmatch.fnmatch(file, 'atom.*.out'):
                self.atoms.append(rhatom.atoms(self.geometry, \
                                               path=('{0}/'+file).format(rhdir)))
        self.Natom = len(self.atoms)

        self.opacity = rhopacity.opacity(self.inputs, self.geometry,\
                                         self.atmos, self.spectrum, \
                                         self.rhdir)
