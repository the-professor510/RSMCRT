import numpy as np
import matplotlib.pyplot as plt
import sys

def plot(radius, count, nPackets, dectType):
    print(count)        
    totalCounts = sum(count)
    print("Detector Type : " + dectType)
    print("Total Diffuse : " + str(totalCounts/nPackets))


    fig = plt.figure(1)
    ax1 = fig.add_subplot()

    ax1.plot(radius, count)
    ax1.set_xlabel("Radius (m)")
    ax1.set_ylabel("Counts (Arb. Units)")
    ax1.set_title(dectType + " Detector")
    plt.show()



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
        
        plot(radius, count, nPackets, "Circle")
            
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
            
        plot(radius, count, nPackets, "Fibre")
            
    case _:     #Default case
        print("Unknown Detector Type")
        sys.exit()
