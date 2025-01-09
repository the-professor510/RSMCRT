import numpy as np
import matplotlib.pyplot as plt
import sys

totalCounts = []
aperture = []
focalLength = 2.0

for j in range(1,11):
    folderName = "RSMCRT/data/detectors/"
    filename = folderName + "detector_" + str(j) + ".dat"

    data = np.fromfile(file=filename, dtype=np.float64, sep="")
    detectorType = data[0]
    match detectorType:
        case 2:     #Detector Type of Fibre
            nPackets = data[1]
            pos = [data[2], data[3], data[4]]
            dir = [data[5], data[6], data[7]]
                        
            focalLength1 = data[8]
            focalLength2 = data[9]
            f1Aperture = data[10]
            f2Aperture = data[11]
            frontOffset = data[12]
            backOffset = data[13]
            frontToPinSep = data[14]
            pinToBackSep = data[15]
            pinAperture = data[16]
            acceptAngle = data[17]
            coreDiameter = data[18]
            
            numBins = (len(data) - 18)/2
            
            radius = []
            count = []
            
            for i in range(19, len(data), 2):
                radius.append(data[i])
                count.append(data[i+1])
                
        case _:     #Default case
            print("Unknown Detector Type")
            sys.exit()
            
    totalCounts.append(sum(count)/nPackets)
    aperture.append(j*0.5)
    print("Total Diffuse " + str(j) + " : " + str(totalCounts[j-1]))


fig = plt.figure(1)
ax1 = fig.add_subplot()

FineAperture = np.linspace(0,5, 1000)
theoreticalTotalCounts = 0.5 * (1- np.cos(np.arctan(FineAperture/focalLength)))
ax1.plot(FineAperture, theoreticalTotalCounts, label = "theoretical")

ax1.scatter(aperture, totalCounts, marker = "x", color = "darkorange", label = "simulated")
ax1.set_xlabel("Radius of lens (m)")
ax1.set_ylabel("Collection Efficiency")
ax1.legend()

plt.show()
