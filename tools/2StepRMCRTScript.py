import subprocess
from argparse import ArgumentParser
import time
import sys
import os
import numpy as np
import read_nrrd_class
import nrrd

if sys.version_info >= (3, 11):
    import tomllib as tomli
else:
    import tomli
# in the future look a pydantic

tic = time.time()
#Run this program from /RSMCRT folder


def gen_command(fpmType, include_Escape, include_Survival, include_Pathlength, toml_File):
    flags = []
    
    #Apply logic to toggle specific flags
    if include_Escape:
        flags.append("-DescapeFunction")

    if include_Survival:
        flags.append("-DsurvivalBias")
    
    if include_Pathlength:
        flags.append("-Dpathlength")

    #Join the flags into a single string for the --flag argument
    flag_string = " ".join(flags)

    #Construct the full command list
    command = [
        "fpm", 
        fpmType, 
        "--flag", flag_string, 
        "--", toml_File
    ]
    return command

#use argument parser to allow the files to be chosen at run time rather than having to save and edit this script each time
parser = ArgumentParser()
parser.add_argument('-f', '--forward_File', help="file path to .toml used to calculate excitation fluence field")
parser.add_argument('-e', '--escape_File', help="file path to .toml used to calculate escape fluence fields for each detector")
parser.add_argument('-fsb', '--forward_SB', default="False", choices=["False", "True"], help="Use survival bias variance reduction when calculating the forwards fluence field")
parser.add_argument('-fpl', '--forward_PL', default="False", choices=["False", "True"], help="Use pathlength variance reduction when calculating the forwards fluence field")
parser.add_argument('-esb', '--escape_SB', default="False", choices=["False", "True"], help="Use survival bias variance reduction when calculating the escape fluence fields")
parser.add_argument('-epl', '--escape_PL', default="False", choices=["False", "True"], help="Use pathlength variance reduction when calculating the escape fluence fields")
parser.add_argument('-ft', '--forward_fpmType', default = "@debugmp", choices=["@debugmp", "@debug", "@runmp"], help="FPM command type for running the forwards code")
parser.add_argument('-at', '--escape_fpmType', default = "@debugmp", choices=["@debugmp", "@debug", "@runmp"], help="FPM command type for running the escape codes")
args = parser.parse_args()

forward_File = "default.toml"
escape_File = "default.toml"
fsb = False
fpl = False
esb = False
epl = False
f_FPM_Type = "@debugmp"
a_FPM_Type = "@debugmp"

#read in arguments
if args.forward_File:
    forward_File = args.forward_File
if args.escape_File:
    escape_File = args.forward_File
if args.forward_SB:
    if args.forward_SB == "True":
        fsb = True
if args.forward_PL:
    if args.forward_PL == "True":
        fpl = True
if args.escape_SB:
    if args.escape_SB == "True":
        esb = True
if args.escape_PL:
    if args.escape_PL == "True":
        epl = True
if args.forward_fpmType:
    f_FPM_Type = args.forward_fpmType
if args.escape_fpmType:
    a_FPM_Type = args.escape_fpmType
    



#ADD READING OF TOML AND TO GET FILE NAMES AND THEN CALCULATION OF THE RAMANDECTEFF

try:
    with open(os.sep.join(["res",forward_File]), "rb") as f:
                fDict = tomli.load(f)
    f.close()
except:
    print(f"Error Failed to Load the forward file: {forward_File}")
    sys.exit(0)

try:
    with open(os.sep.join(["res",escape_File]), "rb") as f:
                eDict = tomli.load(f)
    f.close()
except:
    print(f"Error Failed to Load the escape file: {escape_File}")
    sys.exit(0)
    
fKeys = fDict.keys()
eKeys = eDict.keys()  


#check that the grids are the same, else throw an error
GridMatch = True
fnxg=200
fnyg=200
fnzg=200
fxmax=1.0
fymax=1.0
fzmax=1.0
if ("grid" in fKeys):
    gridKeys = fDict["grid"].keys()
    if ("nxg" in gridKeys):
        fnxg = fDict["grid"]["nxg"]
    if ("nyg" in gridKeys):
        fnyg = fDict["grid"]["nyg"]
    if ("nzg" in gridKeys):
        fnzg = fDict["grid"]["nzg"]
    if ("xmax" in gridKeys):
        fxmax = fDict["grid"]["xmax"]
    if ("ymax" in gridKeys):
        fymax = fDict["grid"]["ymax"]
    if ("zmax" in gridKeys):
        fzmax = fDict["grid"]["zmax"]
    
enxg=200
enyg=200
enzg=200
exmax=1.0
eymax=1.0
ezmax=1.0
if ("grid" in eKeys):
    gridKeys = eDict["grid"].keys()
    if ("nxg" in gridKeys):
        enxg = eDict["grid"]["nxg"]
    if ("nyg" in gridKeys):
        enyg = eDict["grid"]["nyg"]
    if ("nzg" in gridKeys):
        enzg = eDict["grid"]["nzg"]
    if ("xmax" in gridKeys):
        exmax = eDict["grid"]["xmax"]
    if ("ymax" in gridKeys):
        eymax = eDict["grid"]["ymax"]
    if ("zmax" in gridKeys):
        ezmax = eDict["grid"]["zmax"]

if fnxg != enxg:
    print("Error the Grids do not match, check nxg")
    sys.exit(0)
elif fnyg != enyg:
    print("Error the Grids do not match, check nyg")
    sys.exit(0)
elif fnzg != enzg:
    print("Error the Grids do not match, check nzg")
    sys.exit(0)
elif fxmax != exmax:
    print("Error the Grids do not match, check xmax")
    sys.exit(0)
elif fymax != eymax:
    print("Error the Grids do not match, check ymax")
    sys.exit(0)
elif fzmax != ezmax:
    print("Error the Grids do not match, check zmax")
    sys.exit(0)
else:
    print("Grids Checked for Differences")

#check that the GeometryMatches
GeometryMatch = True

fgeom_name = "sphere"
fnumOptProp = 1
fposition = np.array([0.0,0.0,0.0])
fboundingBox = np.array([2.0,2.0,2.0])

if ("geometry" in fKeys):
    geomKeys = fDict["geometry"].keys()
    if "geom_name" in geomKeys:
        fgeom_name = fDict["geometry"]["geom_name"]
    if "numOptProp" in geomKeys:
        fnumOptProp = fDict["geometry"]["numOptProp"]
    if "position" in geomKeys:
        fposition = np.array(fDict["geometry"]["position"])
    if "boundingBox" in geomKeys:
        fboundingBox = np.array(fDict["geometry"]["boundingBox"])
    if fgeom_name == "sphere":
        if "sphereRadius" in geomKeys:
            fsphereRadius = fDict["geometry"]["sphereRadius"]
        else:
            fsphereRadius = 1.0
    if fgeom_name == "box":
        if "BoxDimensions" in geomKeys:
            fBoxDimensions = np.array(fDict["geometry"]["BoxDimensions"])
        else:
            fBoxDimensions = np.array([1.0,1.0,1.0])
    if fgeom_name == "egg":
        if "BottomSphereRadius" in geomKeys:
            fBottomSphereRadius = fDict["geometry"]["BottomSphereRadius"]
        else:
            fBottomSphereRadius = 3.0
        if "TopSphereRadius" in geomKeys:
            fTopSphereRadius = fDict["geometry"]["TopSphereRadius"]
        else:
            fTopSphereRadius = 3.0 * np.sqrt(2.0 - np.sqrt(2.0))
        if "SphereSep" in geomKeys:
            fSphereSep = fDict["geometry"]["SphereSep"]
        else:
            fSphereSep = 3.0 * np.sqrt(2.0 - np.sqrt(2.0))
        if "ShellThickness" in geomKeys:
            fShellThickness = fDict["geometry"]["ShellThickness"]
        else:
            fShellThickness = 0.05
        if "YolkRadius" in geomKeys:
            fYolkRadius =fDict["geometry"]["YolkRadius"]
        else:
            fYolkRadius = 1.5
    if fgeom_name == "cuvette":
        if "outerCuvetteDimensions" in geomKeys:
            fouterCuvetteDimensions = np.array(fDict["geometry"]["outerCuvetteDimensions"])
        else:
            fouterCuvetteDimensions = np.array([1.0,1.0,1.0])
        if "innerCuvetteDimensions" in geomKeys:
            finnerCuvetteDimensions = np.array(fDict["geometry"]["innerCuvetteDimensions"])
        else:
            finnerCuvetteDimensions = np.array([1.0,1.0,1.0])
    if fgeom_name == "multilayerSlab":
        if "thickness" in geomKeys:
            fthickness = fDict["geometry"]["thickness"]
        else:
            fthickness = 1.25*np.ones(fnumOptProp)
        if "xDimensionSize" in geomKeys:
            fxDimensionSize = fDict["geometry"]["xDimensionSize"]
        else:
            fxDimensionSize = 2.0
        if "yDimensionSize" in geomKeys:
            fyDimensionSize = fDict["geometry"]["yDimensionSize"]
        else:
            fyDimensionSize = 2.0

egeom_name = "sphere"
enumOptProp = 1
emur = np.array([1.0])
eposition = np.array([0.0,0.0,0.0])
eboundingBox = np.array([2.0,2.0,2.0])
emur = np.array([1.0])

if ("geometry" in eKeys):
    geomKeys = eDict["geometry"].keys()
    if "geom_name" in geomKeys:
        egeom_name = eDict["geometry"]["geom_name"]
    if "numOptProp" in geomKeys:
        enumOptProp = eDict["geometry"]["numOptProp"]
    if "mur" in geomKeys:
        emur = np.array(eDict["geometry"]["mur"])
    if "position" in geomKeys:
        eposition = np.array(eDict["geometry"]["position"])
    if "boundingBox" in geomKeys:
        eboundingBox = np.array(eDict["geometry"]["boundingBox"])
    if egeom_name == "sphere":
        if "sphereRadius" in geomKeys:
            esphereRadius = eDict["geometry"]["sphereRadius"]
        else:
            esphereRadius = 1.0
    if egeom_name == "box":
        if "BoxDimensions" in geomKeys:
            eBoxDimensions = np.array(eDict["geometry"]["BoxDimensions"])
        else:
            eBoxDimensions = np.array([1.0,1.0,1.0])
    elif egeom_name == "egg":
        if "BottomSphereRadius" in geomKeys:
            eBottomSphereRadius = eDict["geometry"]["BottomSphereRadius"]
        else:
            eBottomSphereRadius = 3.0
        if "TopSphereRadius" in geomKeys:
            eTopSphereRadius = eDict["geometry"]["TopSphereRadius"]
        else:
            eTopSphereRadius = 3.0 * np.sqrt(2.0 - np.sqrt(2.0))
        if "SphereSep" in geomKeys:
            eSphereSep = eDict["geometry"]["SphereSep"]
        else:
            eSphereSep = 3.0 * np.sqrt(2.0 - np.sqrt(2.0))
        if "ShellThickness" in geomKeys:
            eShellThickness = eDict["geometry"]["ShellThickness"]
        else:
            eShellThickness = 0.05
        if "YolkRadius" in geomKeys:
            eYolkRadius =eDict["geometry"]["YolkRadius"]
        else:
            eYolkRadius = 1.5
    elif egeom_name == "cuvette":
        if "outerCuvetteDimensions" in geomKeys:
            eouterCuvetteDimensions = np.array(eDict["geometry"]["outerCuvetteDimensions"])
        else:
            eouterCuvetteDimensions = np.array([1.0,1.0,1.0])
        if "innerCuvetteDimensions" in geomKeys:
            einnerCuvetteDimensions = np.array(eDict["geometry"]["innerCuvetteDimensions"])
        else:
            einnerCuvetteDimensions = np.array([1.0,1.0,1.0])
    elif egeom_name == "multilayerSlab":
        if "thickness" in geomKeys:
            ethickness = eDict["geometry"]["thickness"]
        else:
            ethickness = 1.25*np.ones(enumOptProp)
        if "xDimensionSize" in geomKeys:
            exDimensionSize = eDict["geometry"]["xDimensionSize"]
        else:
            exDimensionSize = 2.0
        if "yDimensionSize" in geomKeys:
            eyDimensionSize = eDict["geometry"]["yDimensionSize"]
        else:
            eyDimensionSize = 2.0
            


if egeom_name != fgeom_name:
    print("Error the Geometries do not match, check geom_name")
    sys.exit(0)
if enumOptProp != fnumOptProp:
    print("Error the Geometries do not match, check numOptProp")
    sys.exit(0)
if (not np.array_equal(eposition, fposition)):
    print("Error the Geometries do not match, check position")
    sys.exit(0)
if (not np.array_equal(eboundingBox, fboundingBox)):
    print("Error the Geometries do not match, check boundingBox")
    sys.exit(0)
if egeom_name == "sphere":
    if esphereRadius != fsphereRadius:
        print("Error the Geometries do not match, check sphereRadius")
        sys.exit(0)
elif egeom_name == "box":
    if (not np.array_equal(eBoxDimensions, fBoxDimensions)):
        print("Error the Geometries do not match, check BoxDimensions")
        sys.exit(0)
elif egeom_name == "egg":
    if eBottomSphereRadius != fBottomSphereRadius:
        print("Error the Geometries do not match, check BottomSphereRadius")
        sys.exit(0)
    if eTopSphereRadius != fTopSphereRadius:
        print("Error the Geometries do not match, check TopSphereRadius")
        sys.exit(0)
    if eSphereSep != fSphereSep:
        print("Error the Geometries do not match, check SphereSep")
        sys.exit(0)
    if eShellThickness != fShellThickness:
        print("Error the Geometries do not match, check ShellThickness")
        sys.exit(0)
    if eYolkRadius != fYolkRadius:
        print("Error the Geometries do not match, check YolkRadius")
        sys.exit(0)
elif egeom_name == "cuvette":
    if (not np.array_equal(eouterCuvetteDimensions, fouterCuvetteDimensions)):
        print("Error the Geometries do not match, check outerCuvetteDimensions")
        sys.exit(0)
    if (not np.array_equal(einnerCuvetteDimensions, finnerCuvetteDimensions)):
        print("Error the Geometries do not match, check innerCuvetteDimensions")
        sys.exit(0)
elif egeom_name == "multilayerSlab":
    if (not np.array_equal(ethickness, fthickness)):
        print("Error the Geometries do not match, check thickness")
        sys.exit(0)
    if exDimensionSize != fxDimensionSize:
        print("Error the Geometries do not match, check xDimensionSize")
        sys.exit(0)
    if eyDimensionSize != fyDimensionSize:
        print("Error the Geometries do not match, check yDimensionSize")
        sys.exit(0)
        
print("Geometries Checked for Differences")

geometryRendered = False
geometryRenderSize = np.array([200,200,200])
geometry_FileName = None
fluence_FileName = None
#read in necessary data to import the illumination field
if "output" in fKeys:
    outputFKeys = fDict["output"].keys()
    if "fluence" in outputFKeys:
        fluence_FileName = os.sep.join(["data", "jmean", fDict["output"]["fluence"]])
    else:
        fluence_FileName = os.sep.join(["data", "jmean", "fluence.nrrd"])
    #end if
        
    if "render_geometry" in outputFKeys:
        geometryRendered = fDict["output"]["render_geometry"]
        if geometryRendered:
            if "render_geometry_name" in outputFKeys:
                geometry_FileName = os.sep.join(["data", fDict["output"]["render_geometry_name"]])
            else:
                geometry_FileName = os.sep.join(["data", "geom_render.nrrd"])
            #end if
            if "render_size" in outputFKeys:
                geometryRenderSize = np.array(fDict["output"])
        #end if
    #end if
#end if

#read in necessary data to import the escape function fields
#first check if geometry has been rendered
if "output" in eKeys:
    outputEKeys = eDict["output"].keys()
    if "render_geometry" in outputEKeys:
        geometryRendered = eDict["output"]["render_geometry"]
        if geometryRendered:
            if "render_geometry_name" in outputFKeys:
                geometry_FileName = os.sep.join(["data", fDict["output"]["render_geometry_name"]])
            else:
                geometry_FileName = os.sep.join(["data", "geom_render.nrrd"])
            #end if
            if "render_size" in outputFKeys:
                geometryRenderSize = np.array(fDict["output"])

if not geometryRendered:
    print(f"WARNING No geometry will be rendered, assume that Raman Conversion Efficiency is Uniform")
elif (not np.array_equal(geometryRenderSize, np.array([enxg, enyg, enzg]))):
    print("WARNING The rendered geometry grid resolution does not match the fluence and escape resoluitons.")
    while True:
            decision = input("Do you want to continue with the assumption that Raman Conversion Efficiency is Uniform? (Y) or (N): ")
            if type(decision) == str:
                if decision.upper().strip() == "Y":
                    geometryRendered = False
                    break
                elif decision.upper().strip() == "N":
                    print("Stopping")
                    sys.exit(0)
            print("Please enter ONLY (Y) or (N)")
     

dectList = None
dectIDList = []
escape_FileName_List = []
if "detectors" in eKeys:
    dectList = eDict["detectors"]
    
    if len(dectList) == 0:
        print("Error No Detectors In Escape TOML Configuration File")
        sys.exit(0)
    
    for index, detector in enumerate(dectList):
        tempKeys = detector.keys()
        if "ID" in tempKeys:
            ID = detector["ID"]
            fileName = os.sep.join(["data", "escape", f"dectID_{ID.strip()}__escape{index+1}.nrrd"])
            dectIDList.append(ID.strip())
            escape_FileName_List.append(fileName)
        else:
            print("Error no ID field. Must specify detector ID for detector with properties:")
            print(detector)
            sys.exit(0)
    
print("Outputs and Detectors checked")  
     



#generate the commands
forward_CMD = gen_command(f_FPM_Type, False, fsb, fpl, forward_File)
escape_CMD = gen_command(a_FPM_Type,  True, esb, epl, escape_File)

print(f"Forward Command: {forward_CMD}")
print(f"Escape Command:  {escape_CMD}")

#"""
#Run the forward code
print("\nRun the Forwards code \n")


try:
    # run() executes the command and waits for it to finish
    result = subprocess.run(
        forward_CMD, 
        capture_output=False,   # Captures stdout and stderr
        text=True,              # Returns output as strings instead of bytes
        check=True              # Raises an error if the Fortran program fails
    )

except subprocess.CalledProcessError as e:
    print(f"Error occurred! Exit code: {e.returncode}")
    print(f"Error message:\n{e.stderr}")
    sys.exit(0)


print("\nRun the Escape code \n")
#Run the escape code
try:
    # run() executes the command and waits for it to finish
    result = subprocess.run(
        escape_CMD, 
        capture_output=False,   # Captures stdout and stderr
        text=True,              # Returns output as strings instead of bytes
        check=True              # Raises an error if the Fortran program fails
    )

except subprocess.CalledProcessError as e:
    print(f"Error occurred! Exit code: {e.returncode}")
    print(f"Error message:\n{e.stderr}")
    sys.exit(0)
print("\nsuccessfully run both simulations\n")
#"""




if geometryRendered:
    plot_nrrd_object = read_nrrd_class.read_nrrd_class()
    gridGeom, hdrGeom = plot_nrrd_object.read_nrrd(geometry_FileName)
    

plot_nrrd_object = read_nrrd_class.read_nrrd_class()
gridFlue, hdrFlue = plot_nrrd_object.read_nrrd(fluence_FileName)


#for calculating the RamanDectEfficiency take into account the relative raman production efficiency between each layer
if (egeom_name == "sphere") or (egeom_name == "box") or (egeom_name == "egg") or (egeom_name == "cuvette") or (egeom_name == "multilayerSlab"):
    layerID = np.arange(1,enumOptProp+1)
    mur = np.append(emur, 0.0)
else:
    layerID = None
    mur = None
    

for index, (escape_FileName, detectorID) in enumerate(zip(escape_FileName_List, dectIDList)):
    plot_nrrd_object = read_nrrd_class.read_nrrd_class()
    gridEscape, hdrEscape = plot_nrrd_object.read_nrrd(escape_FileName)
    
    #write out data as nrrd
    
    
    fileName = f"dectID_{detectorID}__RamanDectEff.nrrd".strip()
    folderRamDectEff = os.sep.join(["data", "RamanDectEff2"])
    write_FileName = os.sep.join([folderRamDectEff, fileName])
    
    
    isExist = os.path.exists(folderRamDectEff)
    isFile = os.path.isfile(write_FileName)
    if not isExist:
        os.makedirs(folderRamDectEff)
    if isFile:
        #the file already exits
        print(f"Warning the file [{fileName}] already exists")
        while True:
            decision = input("Do you want to continue? Enter (Y) or (N): ")
            if type(decision) == str:
                if decision.upper().strip() == "Y":
                    break
                elif decision.upper().strip() == "N":
                    print("Stopping")
                    sys.exit(0)
            
            print("Please enter ONLY (Y) or (N)")
    
    #multiply Fluence by gridEscape to get the Raman Detection Efficiency
    RamanDectEff = np.multiply(np.array(gridFlue) ,np.array(gridEscape))
    
    #Account for different Raman cross sections of different materials
    if geometryRendered:
        for layer, murVal in zip(layerID, mur):
            RamanDectEff[gridGeom==layer] = murVal*RamanDectEff[gridGeom==layer]
    
    #swap the axis
    RamanDectEff = np.reshape(RamanDectEff, (hdrFlue["sizes"][2],hdrFlue["sizes"][1],hdrFlue["sizes"][0]))
                        
    nrrd.write(write_FileName, RamanDectEff, header=hdrFlue, index_order="C")
    print(f"Successfully written out: {write_FileName}")




toc = time.time()
elapsed = int(toc - tic)
hrs, remainder = divmod(elapsed, 3600)
mins, seconds = divmod(remainder, 60)
print(f"\nTotal Time to run: {hrs:02}hrs {mins:02}mins {seconds:.2f}s")

