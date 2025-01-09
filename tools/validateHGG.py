import numpy as np
import matplotlib.pyplot as plt
import sys


folderName = "RSMCRT/data/detectors/"
filename = folderName + "detector_1.dat"

data = np.fromfile(file=filename, dtype=np.float64, sep="")
detectorType = data[0]
match detectorType:
    case 1:     #Detector Type of Cricle
        nPackets = data[1]
        radius = data[2]
        pos = [data[3], data[4], data[5]]
        dir = [data[6], data[7], data[8]]
        numBins = (len(data) - 8)/2
        
        radius = []
        count = []
        
        for i in range(9, len(data), 2):
            radius.append(data[i])
            count.append(data[i+1])
            
    case _:     #Default case
        print("Unknown Detector Type")
        sys.exit()
        
totalCounts = sum(count)/nPackets

# Validating against https://doi.org/10.1016/0169-2607(95)01640-F
# The total diffuse reflection and tranmission for a finite refractive matched slab
print("Theoretical Total Diffure Reflection : 0.09739")
print("Simulated Total Diffuse Reflection : " + str(totalCounts))
print("%Diff : "+ str((np.abs(totalCounts-0.09739)*100)/((0.09739 + totalCounts)/2)))
print()

folderName = "RSMCRT/data/detectors/"
filename = folderName + "detector_2.dat"

data = np.fromfile(file=filename, dtype=np.float64, sep="")
detectorType = data[0]
match detectorType:
    case 1:     #Detector Type of Cricle
        nPackets = data[1]
        radius = data[2]
        pos = [data[3], data[4], data[5]]
        dir = [data[6], data[7], data[8]]
        numBins = (len(data) - 8)/2
        
        radius = []
        count = []
        
        for i in range(9, len(data), 2):
            radius.append(data[i])
            count.append(data[i+1])
            
    case _:     #Default case
        print("Unknown Detector Type")
        sys.exit()
        
totalCounts = sum(count)/nPackets

print("Theoretical Total Diffure Transmission : 0.66096")
print("Simulated Total Diffuse Transmission : " + str(totalCounts))
print("%Diff : "+ str((np.abs(totalCounts-0.66096)*100)/((0.66096 + totalCounts)/2)))