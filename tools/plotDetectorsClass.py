import numpy as np
import matplotlib.pyplot as plt
import sys

class plotDetectorsClass:
       
    def plot(self, radius, count, nPackets, dectType):
        print(radius)
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


    def read_1D_Detector(self, filename):
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
                
                return radius, count, nPackets, numBins, pos, dir, "Circular", \
                    [radius]
                    
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
                    
                
                return radius, count, nPackets, numBins, pos, dir, "Fibre", \
                    [focalLength1, focalLength2, f1Aperture, f2Aperture, \
                        frontOffset, backOffset, frontToPinSep, pinToBackSep, \
                        pinAperture, acceptAngle, coreDiameter] 
                    
            case 3:     #Detector Type of Annulus
                nPackets = data[1]
                radius1 = data[2]
                radius2 = data[3]
                pos = [data[4], data[5], data[6]]
                dir = [data[7], data[8], data[9]]
                numBins = (len(data) - 9)/2
                
                radius = []
                count = []
                
                for i in range(10, len(data), 2):
                    radius.append(data[i])
                    count.append(data[i+1])
                
                return radius, count, nPackets, numBins, pos, dir, "Annulus", \
                    [radius1, radius2]
            
            case _:     #Default case
                print("Unknown Detector Type")
                sys.exit()
