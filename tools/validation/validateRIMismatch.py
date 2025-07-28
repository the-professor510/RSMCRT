#add parent folder to PATH to allow pandalab_base to be imported
import os
import sys

child_dir = os.path.dirname(__file__)
parent_dir = os.path.abspath(os.path.join(child_dir, '..'))
sys.path.append(parent_dir)

import numpy as np
import plotDetectorsClass


folderName = "data/detectors/"
filename = folderName + "detector_1.dat"

readDetectors = plotDetectorsClass.plotDetectorsClass()
radius, count, dectID, nPackets, numBins, pos, dir, typeOfDect, _ = readDetectors.read_1D_Detector(filename)
totalCounts = sum(count)/nPackets

# Validating against https://doi.org/10.1016/0169-2607(95)01640-F
# The total diffuse reflection infinite refractive mismatched slab
print(f"Theoretical Total Diffure Reflection : {0.2600:.5f}")
print(f"Simulated Total Diffuse Reflection : {totalCounts:.5f}")
print(f"%Diff : {((np.abs(totalCounts-0.2600)*100)/((0.2600 + totalCounts)/2)):.5f}")
print()



