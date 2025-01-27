# Hal
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("i")

args = parser.parse_args()

#Given two grids of the same size this will multiply both grids piecewise to find the raman detection efficiency
import numpy as np
import matplotlib.pyplot as plt
import read_nrrd_class
import nrrd
import sys
import os

#Excitation light distribution
#Try either absorb or fluence, comment out which one you don't use
#Personally I only use fluence, if the absorption cross section was set to 0.0

#folderNameAbsorb = "RSMCRT/data/absorb/"
#filenameAbsorb = "absorb.dat"
#file = folderNameAbsorb + filenameAbsorb
#plot_nrrd_object = read_nrrd_class.read_nrrd_class()
#gridAbsorb, hdrAbsorb = plot_nrrd_object.read_nrrd(file)

folderNameFluence = "RSMCRT/data/jmean/"
filenameFluence = "fluence.nrrd"
file = folderNameFluence + filenameFluence
plot_nrrd_object = read_nrrd_class.read_nrrd_class()
gridAbsorb, hdrAbsorb = plot_nrrd_object.read_nrrd(file)


#Escape function for a given detector
folderNameEscape = "RSMCRT/data/escape/"
#i = 11
i = args.i
fileNameEscape = "dectID_Offset"+str(i) + "mm__escape"+str(i)+".nrrd"
file = folderNameEscape + fileNameEscape
plot_nrrd_object = read_nrrd_class.read_nrrd_class()
gridEscape, hdrEscape = plot_nrrd_object.read_nrrd(file)

#write out data as nrrd
folderRamDectEff = "RSMCRT/data/RamanDectEff"
extraID = " "
fileName = f"/{hdrEscape["dector"]}{extraID}.nrrd".strip()
#fileName = "/test.nrrd"
file = folderRamDectEff + fileName
isExist = os.path.exists(folderRamDectEff)
isFile = os.path.isfile(file)
if not isExist:
    os.makedirs(folderRamDectEff)
if isFile:
   #error the file already exits
   print("Error the file already exists")
   sys.exit() 


#check that they are the same size
if(hdrAbsorb["sizes"][0] != hdrEscape["sizes"][0]):
    sys.exit()
elif(hdrAbsorb["sizes"][1] != hdrEscape["sizes"][1]):
    sys.exit(0)
elif(hdrAbsorb["sizes"][2] != hdrEscape["sizes"][2]):
    sys.exit(0)
#else:
    #they are the same size, do nothing

RamanDectEff = np.multiply(np.array(gridAbsorb) ,np.array(gridEscape))

#swap the axis
RamanDectEff = np.reshape(RamanDectEff, (hdrAbsorb["sizes"][2],hdrAbsorb["sizes"][1],hdrAbsorb["sizes"][0]))

nrrd.write(file, RamanDectEff, header=hdrAbsorb, index_order="C")

