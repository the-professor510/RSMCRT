import numpy as np
import matplotlib.pyplot as plt
import plotDetectorsClass

totalCounts = []
aperture = []
focalLength = 2.0

for j in range(1,11):
    folderName = "RSMCRT/data/detectors/"
    filename = folderName + "detector_" + str(j) + ".dat"

    readDetectors = plotDetectorsClass.plotDetectorsClass()
    radius, count, nPackets, numBins, pos, dir, typeOfDect, _ = readDetectors.read_1D_Detector(filename)
                
    totalCounts.append(sum(count)/nPackets)
    aperture.append(j*0.5)
    print(f"Total Diffuse {j} : {(totalCounts[j-1]):.5f}")


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
