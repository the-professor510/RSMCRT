#add parent folder to PATH to allow pandalab_base to be imported
import os
import sys

child_dir = os.path.dirname(__file__)
parent_dir = os.path.abspath(os.path.join(child_dir, '..'))
sys.path.append(parent_dir)

import numpy as np
import matplotlib.pyplot as plt
import read_nrrd_class

folderName = "data/absorb/"
fileName = "absorbValidation2.nrrd"
file = folderName + fileName
plot_nrrd_object = read_nrrd_class.read_nrrd_class()
grid, hdr = plot_nrrd_object.read_nrrd(file)  

folderName2 = "data/absorb/"
fileName2 = "absorbValidation3.nrrd"
file2 = folderName2 + fileName2
plot_nrrd_object = read_nrrd_class.read_nrrd_class()
grid2, hdr = plot_nrrd_object.read_nrrd(file2)  

fig, ax = plt.subplots(2,1, figsize= (5,9))

depths = np.linspace(-2.0, 2.0, int(hdr["sizes"][0]))
TheoryDepths = np.linspace(-2.0,2.0,10000)
ymid = int(hdr["sizes"][1]/2)
zmid = int(hdr["sizes"][2]/2)
data = grid
data2 = grid2

fluence = np.mean(np.mean(data, axis =2), axis = 1)/10
fluence2 = np.mean(np.mean(data2, axis =2), axis = 1)

#""" Used for validate 2
c1 = 5.76
k1 = 1.00
c2 = 1.31
k2 = 10.2
delta = 0.047
norm = 0.115
#"""
fittingFunction = norm * (c1* np.exp((TheoryDepths-1.95)*k1/delta) - c2*np.exp((TheoryDepths-1.95)*k2/delta))
ax[0].plot(TheoryDepths, fittingFunction, label = "Theory - 420 nm", color = "tab:blue")
PercentageError = norm * (c1* np.exp((depths-1.95)*k1/delta) - c2*np.exp((depths-1.95)*k2/delta))
PercentageError = np.abs(fluence - PercentageError)/(0.5*(fluence + PercentageError))*100


ax[1].plot(depths, PercentageError, label = "Percentage Error - 420 nm", color = "darkgray", alpha = 0.5, linestyle = "--")


#""" Used for validate 3
c1 = 6.27
k1 = 1.00
c2 = 1.18
k2 = 14.4
delta = 0.261
norm = 0.0151
#"""
fittingFunction = norm * (c1* np.exp((TheoryDepths-1.95)*k1/delta) - c2*np.exp((TheoryDepths-1.95)*k2/delta))
ax[0].plot(TheoryDepths, fittingFunction, label = "Theory - 630 nm", color = "tab:red")
PercentageError = norm * (c1* np.exp((depths-1.95)*k1/delta) - c2*np.exp((depths-1.95)*k2/delta))
PercentageError = np.abs(fluence2 - PercentageError)/(0.5*(fluence2 + PercentageError))*100
ax[1].plot(depths, PercentageError, label = "Percentage Error - 630 nm", color = "gray", alpha = 0.5, linestyle = "solid")
#ax1.set_title("Refractive Index Mismatch Validation")


ax[0].scatter(depths, fluence, label = "SignedMCRT - 420 nm", marker = "x", color = "tab:pink")
ax[0].scatter(depths, fluence2, label = "SignedMCRT - 630 nm", marker = "x", color = "tab:orange")
ax[0].set_xlabel("Depth, z, (cm)")
ax[0].set_ylabel("Normalised Fluence $\\frac{\Psi(z)}{\Psi_{0}}$ (-)")
ax[0].set_xlim([depths[-14],1.6])
ax[0].set_ylim([0, np.max(fluence)*1.2])

ax[1].set_xlabel("Depth, z, (cm)")
ax[1].set_ylabel("Percentage Error (%)")
ax[1].set_xlim([depths[-14],1.6])
ax[1].set_ylim([0, 12.0])

ax[0].legend(loc = 'upper right')
ax[1].legend(loc = 'upper left')
plt.show()

