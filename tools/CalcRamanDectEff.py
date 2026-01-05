# Hal
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("i")

args = parser.parse_args()

#COMMAND TO RUN THE CALCRAMANDECTEFF MULTIPLE TIMES
#  for i in `seq 1 11`; do python3 RSMCRT/tools/CalcRamanDectEff.py $i ; done

#Given two grids of the same size this will multiply both grids piecewise to find the raman detection efficiency
import numpy as np
import matplotlib.pyplot as plt
import read_nrrd_class
import nrrd
import sys
import os

r"\\wsl.localhost\Ubuntu\home\omguser\MCRT\RamanSMCRT\Results\RamanTest\adjoint\geom_render.nrrd"
#geometry
folderNameGeom = "data/"
fileNameGeom = "geom_render.nrrd"
file = folderNameGeom + fileNameGeom
plot_nrrd_object = read_nrrd_class.read_nrrd_class()
gridGeom, hdrGeom = plot_nrrd_object.read_nrrd(file)

#Excitation light distribution
folderNameFluence = "data/jmean/"
filenameFluence = "fluence.nrrd"
file = folderNameFluence + filenameFluence
plot_nrrd_object = read_nrrd_class.read_nrrd_class()
gridFlue, hdrFlue = plot_nrrd_object.read_nrrd(file)


#Escape function for a given detector
folderNameEscape = "data/escape/"
#i = 11
i = args.i
distance = (2*int(i) - 1)/2
distances = [0.75, 1.25, 1.75, 2.25, 2.75, 3.25, 3.75, 4.25, 4.75]
print(distance)
fileNameEscape = "dectID_"+ "{:.1f}".format(distance) + "mm__escape"+str(int(i))+".nrrd"
file = folderNameEscape + fileNameEscape
plot_nrrd_object = read_nrrd_class.read_nrrd_class()
gridEscape, hdrEscape = plot_nrrd_object.read_nrrd(file)

#write out data as nrrd
folderRamDectEff = "data/RamanDectEff"
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
if(hdrFlue["sizes"][0] != hdrEscape["sizes"][0]):
    sys.exit()
elif(hdrFlue["sizes"][1] != hdrEscape["sizes"][1]):
    sys.exit(0)
elif(hdrFlue["sizes"][2] != hdrEscape["sizes"][2]):
    sys.exit(0)
#else:
    #they are the same size, do nothing

RamanDectEff = np.multiply(np.array(gridFlue) ,np.array(gridEscape))

#account for the different Raman cross sections of different materials
#RamanDectEff[gridGeom==1] = 0
#RamanDectEff[gridGeom==3] = 4*RamanDectEff[gridGeom==3]
#RamanDectEff[gridGeom==5] = 1*RamanDectEff[gridGeom==5]

#swap the axis
RamanDectEff = np.reshape(RamanDectEff, (hdrFlue["sizes"][2],hdrFlue["sizes"][1],hdrFlue["sizes"][0]))

nrrd.write(file, RamanDectEff, header=hdrFlue, index_order="C")

