import numpy as np
import matplotlib.pyplot as plt
import read_nrrd_class

folderName = "RSMCRT/data/absorb/"
fileName = "absorb.nrrd"
file = folderName + fileName
plot_nrrd_object = read_nrrd_class.read_nrrd_class()
grid, hdr = plot_nrrd_object.read_nrrd(file)
    

fig = plt.figure(1)
ax1 = fig.add_subplot()

depths = np.linspace(-2.0, 2.0, int(hdr["sizes"][0]))
ymid = int(hdr["sizes"][1]/2)
zmid = int(hdr["sizes"][2]/2)
data = grid

fluence = np.mean(np.mean(data, axis =2), axis = 1)

ax1.plot(depths, fluence, label = "Simulated")
ax1.set_xlabel("Depth (cm)")
ax1.set_ylabel("Fluence (-)")
ax1.set_xlim([depths[-14],1.6])
ax1.set_ylim([0, np.max(fluence)*1.1])

#""" Used for validate 2
c1 = 5.76
k1 = 1.00
c2 = 1.31
k2 = 10.2
delta = 0.047
norm = 0.115
#"""

""" Used for validate 3
c1 = 6.27
k1 = 1.00
c2 = 1.18
k2 = 14.4
delta = 0.261
norm = 0.0151
#"""

fittingFunction = norm * (c1* np.exp((depths-1.95)*k1/delta) - c2*np.exp((depths-1.95)*k2/delta))
ax1.plot(depths, fittingFunction, label = "Reference")
ax1.set_title("Refractive Index Mismatch Validation")

plt.legend()
plt.show()

