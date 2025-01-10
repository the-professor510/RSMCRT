import numpy as np
import plotDetectorsClass


folderName = "RSMCRT/data/detectors/"
filename = folderName + "detector_1.dat"

readDetectors = plotDetectorsClass.plotDetectorsClass()
radius, count, nPackets, numBins, pos, dir, typeOfDect, _ = readDetectors.read_1D_Detector(filename)
totalCounts = sum(count)/nPackets

# Validating against https://doi.org/10.1016/0169-2607(95)01640-F
# The total diffuse reflection and tranmission for a finite refractive matched slab
print(f"Theoretical Total Diffure Reflection : {0.09739}")
print(f"Simulated Total Diffuse Reflection : {totalCounts:.5f}")
print(f"%Diff : {((np.abs(totalCounts-0.09739)*100)/((0.09739 + totalCounts)/2)):.5f}")
print()

folderName = "RSMCRT/data/detectors/"
filename = folderName + "detector_2.dat"

readDetectors = plotDetectorsClass.plotDetectorsClass()
radius, count, nPackets, numBins, pos, dir, typeOfDect, _ = readDetectors.read_1D_Detector(filename)
totalCounts = sum(count)/nPackets

print(f"Theoretical Total Diffure Transmission : {0.66096}")
print(f"Simulated Total Diffuse Transmission : {(totalCounts):.5f}")
print(f"%Diff : {((np.abs(totalCounts-0.66096)*100)/((0.66096 + totalCounts)/2)):.5f}")