import numpy as np
import plotDetectorsClass


folderName = "RSMCRT/data/detectors/"
filename = folderName + "detector_1.dat"

readDetectors = plotDetectorsClass.plotDetectorsClass()
radius, count, nPackets, numBins, pos, dir, typeOfDect, _ = readDetectors.read_1D_Detector(filename)
totalCounts = sum(count)/nPackets

# Validating against https://doi.org/10.1016/0169-2607(95)01640-F
# The total diffuse reflection and tranmission for a finite refractive matched slab
print("Theoretical Total Diffure Reflection : 0.09739")
print("Simulated Total Diffuse Reflection : " + str(totalCounts))
print("%Diff : "+ str((np.abs(totalCounts-0.09739)*100)/((0.09739 + totalCounts)/2)))
print()

folderName = "RSMCRT/data/detectors/"
filename = folderName + "detector_2.dat"

readDetectors = plotDetectorsClass.plotDetectorsClass()
radius, count, nPackets, numBins, pos, dir, typeOfDect, _ = readDetectors.read_1D_Detector(filename)
totalCounts = sum(count)/nPackets

print("Theoretical Total Diffure Transmission : 0.66096")
print("Simulated Total Diffuse Transmission : " + str(totalCounts))
print("%Diff : "+ str((np.abs(totalCounts-0.66096)*100)/((0.66096 + totalCounts)/2)))