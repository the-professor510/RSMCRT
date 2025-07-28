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

#read in the ramanDectEff data for each spatial offset
for i in range(1,10):
    folderName = "../RSMCRT/data/RamanDectEff/"
    fileName = f"{(i/2+0.25):.2f}mm .nrrd"
    file = folderName + fileName
    plot_nrrd_object = read_nrrd_class.read_nrrd_class()
    grid, hdr = plot_nrrd_object.read_nrrd(file)
    
    
    NormalContribution.append(np.sum(grid[gridGeom==2]))
    TumorContribution.append(np.sum(grid[gridGeom==3]))
    spatialOffset.append((i/2+0.25))
    
NormalContribution = np.array(NormalContribution)
TumorContribution = np.array(TumorContribution)
RelTumContribution = TumorContribution/(TumorContribution+NormalContribution)

fig = plt.figure(1)
ax1 = fig.add_subplot()
ax1.plot(spatialOffset, RelTumContribution)
ax1.set_xlim(0,5.0)
ax1.set_ylim(0,0.8)
plt.show()
    
    
    
    
    