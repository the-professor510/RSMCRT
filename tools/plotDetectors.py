import numpy as np
import matplotlib.pyplot as plt
import sys


folderName = "RSMCRT/data/detectors/"
filename = folderName + "detector_2.dat"
NoOpticalDepths = 5
OpticalDepths = [0.1, 1, 10, 30, 100] 

NoTimeBins = 1000
TimeBinResolution = 0.000000001
MaxTime = NoTimeBins*TimeBinResolution
times = np.linspace(0,MaxTime,NoTimeBins, endpoint=False)


data = np.fromfile(file=filename, dtype=np.float64, sep="")
detectorType = data[0]
match detectorType:
    case 1:     #Detector Type of Cricle
        radius = data[1]
        pos = [data[2], data[3], data[4]]
        dir = [data[5], data[6], data[7]]
        numBins = len(data) - 8
        
        radius = []
        count = []
        
        for i in range(8, len(data), 2):
            radius.append(data[i])
            count.append(data[i+1])
            
    case _:     #Default case
        print("Unknown Detector Type")
        sys.exit()
        
#print(numBins)        
#print(radius)
#print(count)
totalCounts = sum(count)

print(count)
#print("Total counts: " + str(totalCounts))
print("Total Diffuse ... : " + str(totalCounts/1000000))


fig = plt.figure(1)
ax1 = fig.add_subplot()

ax1.plot(radius, count)
ax1.set_xlabel("Radius (m)")
ax1.set_ylabel("Counts (Arb. Units)")
ax1.set_title("Radius Detector")
plt.show()