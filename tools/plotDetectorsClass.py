import numpy as np
import matplotlib.pyplot as plt
import sys

class plotDetectorsClass:
       
    def plot(self, radius, count, nPackets, dectType):
             
        totalCounts = sum(count)
        print(f"Detector Type : {dectType}")
        print(f"Total Diffuse : {totalCounts/nPackets}")


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
        n = 1
        match detectorType:
            case 1:     #Detector Type of Cricle
                dectIDLen = data[n]
                dectID = ""
                n+=1
                for i in range(n, n+int(dectIDLen)):
                    dectID = dectID + chr(int(data[i]))
                    n += 1
                nPackets = data[n]
                radius = data[n+1]
                pos = [data[n+2], data[n+3], data[n+4]]
                dir = [data[n+5], data[n+6], data[n+7]]
                numBins = (len(data) - (n+8))/2
                
                radius = []
                count = []
                
                for i in range(n+8, len(data), 2):
                    radius.append(data[i])
                    count.append(data[i+1])
                
                return radius, count, dectID, nPackets, numBins, pos, dir, "Circular", \
                    [radius]
                    
            case 2:     #Detector Type of Fibre
                dectIDLen = data[n]
                dectID = ""
                n+=1
                for i in range(n, n+int(dectIDLen)):
                    dectID = dectID + chr(int(data[i]))
                    n += 1
                nPackets = data[n]
                pos = [data[n+1], data[n+2], data[n+3]]
                dir = [data[n+4], data[n+5], data[n+6]]
                            
                focalLength1 = data[n+7]
                focalLength2 = data[n+8]
                f1Aperture = data[n+9]
                f2Aperture = data[n+10]
                frontOffset = data[n+11]
                backOffset = data[n+12]
                frontToPinSep = data[n+13]
                pinToBackSep = data[n+14]
                pinAperture = data[n+15]
                acceptAngle = data[n+16]
                coreDiameter = data[n+17]
                
                numBins = (len(data) - (n+18))/2
                
                radius = []
                count = []
                
                for i in range(n+18, len(data), 2):
                    radius.append(data[i])
                    count.append(data[i+1])
                    
                
                return radius, count, dectID, nPackets, numBins, pos, dir, "Fibre", \
                    [focalLength1, focalLength2, f1Aperture, f2Aperture, \
                        frontOffset, backOffset, frontToPinSep, pinToBackSep, \
                        pinAperture, acceptAngle, coreDiameter] 
                    
            case 3:     #Detector Type of Annulus
                dectIDLen = data[n]
                dectID = ""
                n+=1
                for i in range(n, n+int(dectIDLen)):
                    dectID = dectID + chr(int(data[i]))
                    n += 1
                nPackets = data[n]
                radius1 = data[n+1]
                radius2 = data[n+2]
                pos = [data[n+3], data[n+4], data[n+5]]
                dir = [data[n+6], data[n+7], data[n+8]]
                numBins = (len(data) - (n+9))/2
                
                radius = []
                count = []
                
                for i in range(n+9, len(data), 2):
                    radius.append(data[i])
                    count.append(data[i+1])
                
                return radius, count, dectID, nPackets, numBins, pos, dir, "Annulus", \
                    [radius1, radius2]
            
            case _:     #Default case
                print("Unknown Detector Type")
                sys.exit()
