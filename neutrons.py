import sys
import sdf
import openpmd_api as api
import time
import math
import numpy as np
from numpy.random import random

SCALAR = api.Mesh_Record_Component.SCALAR
Unit_Dimension = api.Unit_Dimension

fname = "Data/" + sys.argv[1]
varnames = ["Grid_Particles_" + sys.argv[2], "Particles_Px_" + sys.argv[2], "Particles_Py_" + sys.argv[2],
            "Particles_Weight_" + sys.argv[2]]

try:
    d = sdf.read(fname)
except:
    print('File "%s" not found' % fname)
    sys.exit()

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

coords = d.__dict__[varnames[0]]
px = d.__dict__[varnames[1]]
py = d.__dict__[varnames[2]]
wd = d.__dict__[varnames[3]]

w1 = wd.data[0]  # all ions have the same weight

mi = 2 * 1836 * 9.1e-31
elch = 1.6e-19

energy = np.divide(np.add(np.square(px.data), np.square(py.data)), 2 * mi * elch)
# distance from source (foil)
dist = 0.02

de = 10000  # resolution of the spectrum in eV
nb = math.floor(np.max(energy) / de)
counts, binedges = np.histogram(energy, bins=nb)

# this will be needed in the next version
# vels=np.sqrt(np.divide(binedges,mi/2))
# drifttime=dist/vels[0]
# drift=np.substract(np.multiply(vels, drifttime), dist)
# --------------------------------------------------------


# parameters of the target, which is bombarded by the D+ ions
dens = 6.e28  # density in 1/m^3
lt = 10.e-3  # length in m

# cross section of neutron production for 0.5 MeV ion energy:
sigma = np.multiply(np.power(binedges/1000 - 100, 0.3333), 1.0e-30)  # in m^2

rb = 10e-6  # radius of the ion beam
Nn=0
vx=[0]
nw=10000
for i in range(len(sigma)-1):
    if sigma[i]>0:
        px = np.sqrt(2 * mi * binedges[i] * elch)  # momentum (in the x direction) in kg m/s
        Nd = counts[i]*w1*rb**2/nw
        inc = int(round(Nd * dens * lt * sigma[i]))
        if inc>0:
            vx = np.append(vx,np.multiply(np.ones(inc),px/mi))
            Nn = Nn + inc

vx=vx[1:]
print("Number of neutron macroparticles:", Nn)

# Neutron energy: 2.45 MeV
vn = math.sqrt(2 * 2.45e6 * 1.6e-19 / 1.67e-27)

posx = 1.e6 * lt * random(Nn)
r = 1.e6 * rb * random(Nn)  # positions will be saved in units of micron
a = 2 * math.pi * random(Nn)

posy = np.multiply(r, np.cos(a))
posz = np.multiply(r, np.sin(a))

aa = math.pi * random(Nn)
b = 2 * math.pi * random(Nn)

velx = np.add(vn * np.cos(aa), vx)
velr = vn * np.sin(aa)
vely = np.multiply(velr, np.cos(b))
velz = np.multiply(velr, np.sin(b))

id = np.arange(1, Nn + 1)

weight = nw*np.ones(Nn)

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
neutrons.set_attribute("numParticles", Nn)

d = api.Dataset(id.dtype, id.shape)
neutrons["id"][SCALAR].reset_dataset(d)
neutrons["id"][SCALAR].store_chunk(id)

d = api.Dataset(weight.dtype, weight.shape)
neutrons["weight"][SCALAR].reset_dataset(d)
neutrons["weight"][SCALAR].store_chunk(weight)

d = api.Dataset(posx.dtype, posx.shape)
neutrons["position"]["x"].reset_dataset(d)
neutrons["position"]["y"].reset_dataset(d)
neutrons["position"]["z"].reset_dataset(d)
neutrons["position"]["x"].set_unit_SI(1.e-6)
neutrons["position"]["y"].set_unit_SI(1.e-6)
neutrons["position"]["z"].set_unit_SI(1.e-6)
neutrons["position"].set_unit_dimension({Unit_Dimension.L: 1})
neutrons["position"]["x"].store_chunk(posx)
neutrons["position"]["y"].store_chunk(posy)
neutrons["position"]["z"].store_chunk(posz)

d = api.Dataset(velx.dtype, velx.shape)
neutrons["velocity"]["x"].reset_dataset(d)
neutrons["velocity"]["y"].reset_dataset(d)
neutrons["velocity"]["z"].reset_dataset(d)
neutrons["velocity"]["x"].set_unit_SI(1)
neutrons["velocity"]["y"].set_unit_SI(1)
neutrons["velocity"]["z"].set_unit_SI(1)
neutrons["velocity"].set_unit_dimension({Unit_Dimension.L: 1, Unit_Dimension.T: -1})
neutrons["velocity"]["x"].store_chunk(velx)
neutrons["velocity"]["y"].store_chunk(vely)
neutrons["velocity"]["z"].store_chunk(velz)

series.flush()
del series
