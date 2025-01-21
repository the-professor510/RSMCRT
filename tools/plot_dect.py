import numpy as np
import plotDetectorsClass


folderName = "RSMCRT/data/detectors/"
for i in range(1, 12):
    filename = folderName + "detector_" + str(i) + ".dat"

    readDetectors = plotDetectorsClass.plotDetectorsClass()
    radius, count, dectID, nPackets, numBins, pos, dir, typeOfDect, _ = readDetectors.read_1D_Detector(filename)

    print(dectID)
    readDetectors.plot(radius, count, nPackets, typeOfDect)