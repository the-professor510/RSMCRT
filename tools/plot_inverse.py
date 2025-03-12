import numpy as np
import plotInverseClass

folderName = "RSMCRT/data/inverse/"
filename = folderName + "inverse.dat"

readInverse = plotInverseClass.plotInverseClass()
mus, mua, hgg, n, error, bestGuessIndx, bestMus, bestMua, besthgg, bestn, bestError = readInverse.read_1D_Detector(filename)


readInverse.plot1D(mus, error, bestGuessIndx, "mus")
readInverse.plot1D(mua, error, bestGuessIndx, "mua")
readInverse.plot2D(mus, mua, error, bestGuessIndx)

readInverse.plot1D(hgg, error, bestGuessIndx, "hgg")
readInverse.plot1D(n, error, bestGuessIndx, "n")


