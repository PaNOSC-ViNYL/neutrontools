from IonDistribution import IonDist
import matplotlib.pyplot as plt
import numpy as np
from NeutronDistribution import NeutronSource

mydist = IonDist('Data', 10, 'proton')

mydist.getIons()


mydist.genSpectrum(10000)
fig = plt.figure()
ax = plt.axes()

ax.plot(mydist.binedges/1000, np.log10(mydist.counts))
plt.xlabel("Ion Energy (keV)")
plt.ylabel("log10(Number/10 keV)")
plt.show()

neutrons = NeutronSource(0.01, 1.e28, [mydist.binedges, mydist.counts])
neutrons.genNeutrons(10000)

fig1 = plt.figure()
ax1 = plt.axes()

ax1.scatter(neutrons.data[4]/3.e8, neutrons.data[5]/3.e8)
plt.xlabel("vy/mc")
plt.ylabel("vz/mc")
plt.show()

neutrons.saveData()
