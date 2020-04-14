import openpmd_api as api
import time
import datetime
import math
import numpy as np
from numpy.random import random
SCALAR = api.Mesh_Record_Component.SCALAR
Unit_Dimension = api.Unit_Dimension

en= 5.e5    # 0.5 MeV deuterium ions
            # this is an input parameter here, the ion data should be
            # obtained by using the getIons(dat) function, will be defined below
mi=2*1836*9.1e-31
elch=1.6e-19
px=np.sqrt(2*mi*en*elch)            #momentum (in the x direction) in kg m/s
vx=px/mi

Nd=1.e7         #number of ions
rb=10e-6        #radius of the ion beam

#parameters of the target, which is bombarded by the D+ ions

dens=1.e28   #density in 1/m^3
lt=1.e-3     #length in m

#cross section of neutron production for 0.5 MeV ion energy:
sigma=2.e-29   # in m^2

#Number of neutrons generated:

Nn=math.floor(Nd*dens*lt*sigma)  #in this example it is Nn=2000

#Since the target density is uniform, the neutrons will be generated
#uniformly. Their position will be chosen randomly within 0 and lt

#Neutron energy: 2.45 MeV

vn=math.sqrt(2*2.45e6*1.6e-19/1.67e-27)

posx=1.e6*lt*random(Nn)
r=1.e6*rb*random(Nn)      #positions will be saved in units of micron
a=2*math.pi*random(Nn)

posy=np.multiply(r,np.cos(a))
posz=np.multiply(r,np.sin(a))

aa=math.pi*random(Nn)
b=2*math.pi*random(Nn)

velx=(vn+vx)*np.cos(aa)
velr=vn*np.sin(aa)
vely=np.multiply(velr, np.cos(b))
velz=np.multiply(velr, np.sin(b))

id =np.arange(1,Nn+1)

series = api.Series(
    "dataMD.h5",
    api.Access_Type.create)
dateNow = time.strftime('%Y-%m-%d %H:%M:%S %z', time.localtime())
print("Default settings:")
print("basePath: ", series.base_path)
print("openPMD version: ", series.openPMD)
print("iteration format: ", series.iteration_format)

series.set_openPMD("1.1.0")
series.set_openPMD_extension(0)
series.set_author("Zsolt Lecz<zsolt.lecz@eli-alps.hu>")
series.set_particles_path("particles")
series.set_date(dateNow)
series.set_iteration_encoding(api.Iteration_Encoding.group_based)
series.set_software("EPOCH")
series.set_software_version("4.8.3")
#series.set_attribute("forceField","eam/alloy")
#series.set_attribute("forceFieldParameter","pair_coeff * * Cu_mishin1.eam.alloy Cu")

curStep = series.iterations[0]
curStep.set_time(0.0)        .set_time_unit_SI(1e-15)
curStep.set_attribute("step",np.uint64(0))
curStep.set_attribute("stepOffset",np.uint64(0))
curStep.set_attribute("timeOffset",np.float32(0))

neutrons = curStep.particles["neutrons"]

d = api.Dataset(id.dtype, id.shape)
neutrons["id"][SCALAR].reset_dataset(d)
neutrons["id"][SCALAR].store_chunk(id)

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