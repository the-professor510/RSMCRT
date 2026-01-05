#add parent folder to PATH to allow pandalab_base to be imported
import os
import sys

child_dir = os.path.dirname(__file__)
parent_dir = os.path.abspath(os.path.join(child_dir, '..'))
sys.path.append(parent_dir)

import numpy as np
import matplotlib.pyplot as plt
import read_nrrd_class

spatialOffset = []
NormalContribution = []
TumorContribution = []


#geometry
folderNameGeom = "../RSMCRT/data/"
fileNameGeom = "geom_render.nrrd"
file = folderNameGeom + fileNameGeom
plot_nrrd_object = read_nrrd_class.read_nrrd_class()
gridGeom, hdrGeom = plot_nrrd_object.read_nrrd(file)

trueDataX = np.array([0.750, 1.246, 1.743, 2.239, 2.736, 3.242, 3.744, 4.245, 4.747])
trueDataY = np.array([0.260, 0.335, 0.383, 0.460, 0.535, 0.512, 0.562, 0.607, 0.608])
TopError = np.array([0.276, 0.362, 0.425, 0.517, 0.596, 0.606, 0.647, 0.685, 0.646])

#read in the ramanDectEff data for each spatial offset
distances = [0.75, 1.25, 1.75, 2.25, 2.75, 3.25, 3.75, 4.25, 4.75]
for i in range(1,10):
    folderName = "../RSMCRT/data/RamanDectEff/"
    fileName = f"{distances[i-1]:.2f}mm .nrrd"
    file = folderName + fileName
    
    folderName = "../RSMCRT/data/raman/"
    fileName = f"dectID_{distances[i-1]:.2f}mm__escape{i}.nrrd"
    file = folderName + fileName
    
    print(file)
    
    
    plot_nrrd_object = read_nrrd_class.read_nrrd_class()
    grid, hdr = plot_nrrd_object.read_nrrd(file)
    
    
    NormalContribution.append(np.sum(grid[gridGeom==3]))
    TumorContribution.append(np.sum(grid[gridGeom==5]))
    spatialOffset.append(distances[i-1])
    
NormalContribution = np.array(NormalContribution)
TumorContribution = np.array(TumorContribution)
RelTumContribution = TumorContribution/(TumorContribution+NormalContribution)

fig = plt.figure(1)
ax1 = fig.add_subplot()
ax1.plot(spatialOffset, RelTumContribution)
ax1.errorbar(trueDataX, trueDataY, yerr = TopError-trueDataY, marker = "x")
ax1.set_xlim(0,5.0)
ax1.set_ylim(0,0.8)
plt.show()
    
    
    
    
    