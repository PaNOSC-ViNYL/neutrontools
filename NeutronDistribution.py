import sys
import numpy as np
import math
from numpy.random import random
import openpmd_api as api
import time


class NeutronSource:
    __s = []    #cross section calculted for the current ion energy spectrum
    Nn = 0
    __dims = 8
    data = []

    def __init__(self, length, density, array):
        self.lt = length
        self.den = density
        self.spec = array

        if len(array) != 2:
            print("The spectrum should contain both axes vaues.")
            sys.exit()
        self.__s = [1.0e-30 * (en/1000 - 100) ** (1.0 / 3) for en in self.spec[0] if en < 2.e6]

    def genNeutrons(self, nw):
        rb = 10e-6  # radius of the ion beam
        mi = 2 * 1836 * 9.1e-31
        elch = 1.6e-19
        vx = [0]

        for i in range(len(self.__s) - 1):
            if self.__s[i] > 0 and self.spec[0][i]<2.e6:
                px = np.sqrt(2 * mi * self.spec[0][i] * elch)  # momentum (in the x direction) in kg m/s
                Nd = self.spec[1][i] * rb ** 2 / nw
                inc = int(round(Nd * self.den * self.lt * self.__s[i]))
                if inc > 0:
                    vx = np.append(vx, np.multiply(np.ones(inc), px / mi))
                    self.Nn += inc

        vx = vx[1:]
        print("Number of neutron macroparticles:", self.Nn)

        self.data = np.zeros(shape=(self.__dims, self.Nn))
        # Neutron energy: 2.45 MeV
        vn = math.sqrt(2 * 2.45e6 * elch / 1.67e-27)
        self.data[0] = 1.e6 * self.lt * random(self.Nn)
        r = 1.e6 * rb * random(self.Nn)  # positions will be saved in units of micron
        a = 2 * math.pi * random(self.Nn)

        self.data[1] = np.multiply(r, np.cos(a))
        self.data[2] = np.multiply(r, np.sin(a))

        aa = math.pi * random(self.Nn)
        b = 2 * math.pi * random(self.Nn)

        self.data[3] = np.add(vn * np.cos(aa), vx)
        velr = vn * np.sin(aa)
        self.data[4] = np.multiply(velr, np.cos(b))
        self.data[5] = np.multiply(velr, np.sin(b))
        self.data[6] = np.arange(1, self.Nn + 1)   # id
        self.data[7] = nw * np.ones(self.Nn)       # weight

    def saveData(self):
        SCALAR = api.Mesh_Record_Component.SCALAR
        Unit_Dimension = api.Unit_Dimension

        series = api.Series(
            "NeutronData.h5",
            api.Access_Type.create)
        dateNow = time.strftime('%Y-%m-%d %H:%M:%S %z', time.localtime())
        print("Default settings:")
        print("basePath: ", series.base_path)
        print("openPMD version: ", series.openPMD)
        print("iteration format: ", series.iteration_format)

        series.set_openPMD("1.1.0")
        # series.set_openPMD_extension("BeamPhysics;SpeciesType")
        series.set_attribute("openPMDextension", "BeamPhysics;SpeciesType")
        series.set_author("Zsolt Lecz<zsolt.lecz@eli-alps.hu>")
        series.set_particles_path("particles")
        series.set_date(dateNow)
        series.set_iteration_encoding(api.Iteration_Encoding.group_based)
        series.set_software("EPOCH", "4.8.3")
        # series.set_software_version("4.8.3")
        # series.set_attribute("forceField","eam/alloy")
        # series.set_attribute("forceFieldParameter","pair_coeff * * Cu_mishin1.eam.alloy Cu")

        curStep = series.iterations[0]
        curStep.set_time(0.0).set_time_unit_SI(1e-15)
        curStep.set_attribute("step", np.uint64(0))
        curStep.set_attribute("stepOffset", np.uint64(0))
        curStep.set_attribute("timeOffset", np.float32(0))

        neutrons = curStep.particles["neutrons"]
        neutrons.set_attribute("speciesType", "neutron")
        neutrons.set_attribute("numParticles", self.Nn)

        d = api.Dataset(self.data[6].dtype, self.data[6].shape)
        neutrons["id"][SCALAR].reset_dataset(d)
        neutrons["id"][SCALAR].store_chunk(self.data[6])

        d = api.Dataset(self.data[7].dtype, self.data[7].shape)
        neutrons["weight"][SCALAR].reset_dataset(d)
        neutrons["weight"][SCALAR].store_chunk(self.data[7])

        d = api.Dataset(self.data[0].dtype, self.data[0].shape)
        neutrons["position"]["x"].reset_dataset(d)
        neutrons["position"]["y"].reset_dataset(d)
        neutrons["position"]["z"].reset_dataset(d)
        neutrons["position"]["x"].set_unit_SI(1.e-6)
        neutrons["position"]["y"].set_unit_SI(1.e-6)
        neutrons["position"]["z"].set_unit_SI(1.e-6)
        neutrons["position"].set_unit_dimension({Unit_Dimension.L: 1})
        neutrons["position"]["x"].store_chunk(self.data[0])
        neutrons["position"]["y"].store_chunk(self.data[1])
        neutrons["position"]["z"].store_chunk(self.data[2])

        d = api.Dataset(self.data[0].dtype, self.data[0].shape)
        neutrons["velocity"]["x"].reset_dataset(d)
        neutrons["velocity"]["y"].reset_dataset(d)
        neutrons["velocity"]["z"].reset_dataset(d)
        neutrons["velocity"]["x"].set_unit_SI(1)
        neutrons["velocity"]["y"].set_unit_SI(1)
        neutrons["velocity"]["z"].set_unit_SI(1)
        neutrons["velocity"].set_unit_dimension({Unit_Dimension.L: 1, Unit_Dimension.T: -1})
        neutrons["velocity"]["x"].store_chunk(self.data[3])
        neutrons["velocity"]["y"].store_chunk(self.data[4])
        neutrons["velocity"]["z"].store_chunk(self.data[5])

        series.flush()
        del series
