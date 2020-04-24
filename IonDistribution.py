import sys
import sdf
import os
import math
import numpy as np


class IonDist:
    path = ''
    momentum = []
    weight = 0
    counts = []
    binedges = []

    def __init__(self, folder, file_num, name):
        if file_num >= 0:
            self.__file_num = file_num
        else:
            raise ValueError("The parameter 'file_num' must be a positive integer.")
        if isinstance(folder, str):
            if not os.path.isdir(folder):
                raise IOError("The parameters %s is not a valid directory name." % folder)
        else:
            raise TypeError("The first parameter should be a srting.")

        self.pname = name
        if file_num < 10:
            self.path = folder + '/000' + str(file_num) + ".sdf"
        elif file_num < 100:
            self.path = folder + '/00' + str(file_num) + ".sdf"
        elif file_num < 1000:
            self.path = folder + '/0' + str(file_num) + ".sdf"
        else:
            self.path = folder + str(file_num) + ".sdf"

    def getIons(self):
        varnames = ["Grid_Particles_" + self.pname, "Particles_Px_" + self.pname, "Particles_Weight_" + self.pname]
        try:
            d = sdf.read(self.path)
        except:
            print('File "%s" not found' % self.path)
            sys.exit()
        i = 0
        for i in range(len(varnames)):
            if not hasattr(d, varnames[i]):
                print('Variable "%s" not found in file' % varnames[i])
                break

        if i != (len(varnames) - 1):
            if i > 1:
                varnames = varnames[0:i - 1]
            elif i == 1:
                print('Particle data is incomplete: no momentum!')
                sys.exit()
            else:
                print('No particle data!')
                sys.exit()

        #coords = d.__dict__[varnames[0]]
        px = d.__dict__[varnames[1]]
        #py = d.__dict__[varnames[2]]  ... here more components can be included in higher dimensions
        wd = d.__dict__[varnames[-1]]
        self.momentum = px.data
        self.weight = wd.data[1]

    def genSpectrum(self, de):
        mi = 2 * 1836 * 9.1e-31
        elch = 1.6e-19
        if len(self.momentum)==0:
            print("Read the particle data first!")
            sys.exit()
        energy = np.divide(np.square(self.momentum), 2 * mi * elch)
        de = 10000  # resolution of the spectrum in eV
        nb = math.floor(np.max(energy) / de)
        c1, b1 = np.histogram(energy, bins=nb)
        self.binedges = b1[1:]
        self.counts=np.multiply(c1, self.weight)
        print("Number of energy bins: %s" % len(self.binedges))
