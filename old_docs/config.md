---
title: config
---
# Config file settings

The configuration file format used is Tom's Obvious Minimal Language ([TOML](https://toml.io/en/)).
The below sections describe the tables (dictionaries) that are able to be defined for SignedMCRT.

## Source

This table defines the parameters for the light source used in the simulation it can have the following:

| Parameter | Type | Options | Default | Notes |
|:---------:|:----:|:-------:|:-------:|:-----:|
| name | string | point, circular, uniform, pencil, annulus, focus | point | - |
| nphotons | integer | - | 1000000 | - |
| position | float array size 3 | - | [0.0, 0.0, 0.0] | Default value only set for point, annular, and focus source type |
| direction | float array size 3 or string | - | -z | String type applies to circular, uniform and pencil sources |
| point1 | float array size 3 | - | [-1.0, -1.0, -1.0] | Used by uniform source only to set location and size of source |
| point2 | float array size 3 | - | [2.0, 0.0, 0.0] | See Above |
| point3 | float array size 3 | - | [0.0, 2.0, 0.0] | See Above |
| radius | float | - | 0.5 | Used by circular source  |
| rhi | float | - | 0.6 | Annular source upper radius |
| rlo | float | - | 0.5 | Annular source lower radius |
| annulus_type | string | gaussian, tophat | gaussian | Type of annular beam |
| focalLength | float | - | 1.0 | Used by annular and focus, the distance from zmax to the focus point of the annular and focus sources, positive for converging beam, negative for diverging beam |
| rotation | float array size 3 | - | [1.0, 0.0, 0.0] | Beam direction for annular and focus |
| focus_type | string | gaussian, circle, square | gaussian | Shape of focus source |
| beam_size | float | - | 0.5 | size of focus beam, gaussian: 1/e radius, circle: radius, square: half height/width |
| spectrum_type | string | constant, 1D, 2D | constant | Type of spectrum used |
| spectrum_file | string | - | - | filename of 1D or 2D spectrum/image |
| cell_size | float array size 2 | - | - | size of pixel in 2D spectrum in simulation units. |
| wavelength | float | - | 500 nm | Wavelength for constant spectrum |

**Note** point1, point2, and point3 define a rectangle. Point1 is the origin,point2 and point3 are the vectors that describe the sides.

## Grid

| Parameter | Type | Default | Notes |
|:---------:|:----:|:-------:|:-----:|
| nxg | integer | 200 | Number of voxel in x direction |
| nyg | integer | 200 | Number of voxel in y direction |
| nyg | integer | 200 | Number of voxel in z direction |
| xmax | float | 1.0 | Half size of simulated medium in x direction |
| ymax | float | 1.0 | Half size of simulated medium in y direction |
| zmax | float | 1.0 | Half size of simulated medium in z direction |
| units | string | cm | Units of simulation (currently need to manually adjust optical properties to account) |

## Geometry

| Parameter | Type | Options | Default | Notes |
|:---------:|:----:|:-------:|:-------:|:-----:|
| geom_name | string | sphere, box, egg, exp | sphere | Name of experiment for metadata |
| numOptProp | integer | - | 1 | Size of optical property arrays, for egg and box scene must be of size 1, for egg scene must be of size 3 |
| mua | float array size numOptProp | - | 0.0 | Absorption Coefficient |
| mus | float array size numOptProp | - | 1.0 | Scattering Coefficient |
| mur | float array size numOptProp | - | 0.0 | Raman Coefficient |
| hgg | float array size numOptProp | - | 0.0 | Henyey-Greenstein Coefficient |
| n | float array size numOptProp | - | 0.0 | Refractive Index | 
| position | float array size 3 | - | [0.0, 0.0, 0.0] | Position of the sphere, box, and egg | 
| boundingBox | float array size 3 | - | [2.0, 2.0, 2.0] | Dimensions of bounding box with optical properties of a vacuum | 
| sphereRadius | float | - | 1.0 | Radius of Sphere | 
| BoxDimensions | float array size 3 | - | [1.0, 1.0, 1.0] | Dimensions of box | 
| BottomSphereRadius | float | - | 3.0 | Radius of Bottom Sphere of Moss egg |
| SphereSep | float | - | 3.0 * sqrt( 2.0 - sqrt(2.0))  | Seperation between top and bottom sphere of Moss egg, must be equal or less than above |
| TopSphereRadius | float | - | 3.0 * sqrt( 2.0 - sqrt(2.0)) | Radius of Top Sphere of Moss egg, must be equal or less than above |
| ShellThickness | float | 1.0 <-> 0.0 | 0.05 | Corresponds to thickness of the shell, 1.0 corresponds to no shell, 0.0 corresponds to no albumen |
| YolkRadius | float | - | 1.5 | Radius of yolk inside the egg |
| tau | float | - | 10.0 | Tau value for MCRT scattering test experiment |
| num_spheres | integer | - | 10 | Number of random spheres for sphere scene, not supported |
| musb | float | - | 0.00 | Optical properties for experimental geometry for whiskey Raman sensing paper |
| muab | float | - | 0.01 | See above |
| musc | float | - | 0.00 | See Above |
| muac | float | - | 0.01 | See Above |
| hgga | float | - | 0.70 | See Above |


## Detectors

| Parameter | Type | Options | Default | Notes |
|:---------:|:----:|:-------:|:-------:|:-----:|
| type | string | annulus, circle, fibre, camera | - | - |
| position | float array size 3 | - | NO DEFAULT! | Central position of detector |
| direction | float array size 3 | - | [0.0, 0.0, -1.0] | Propagation direction of accepted rays |
| radius | float | - | 1.0 | Radius of circular detector |
| radius1 | float | - | - | Inner radius of annular detector |
| radius2 | float | - | - | Outer radius of annulus detector. Must be larger than radius1 |
| focalLength1 | float | - | 1.0 | Front lens focal length in a 4f system fibre collection system |
| focalLength2 | float | - | 1.0 | Back lens focal length in a 4f system fibre collection system |
| f1Aperture | float | - | 1.0 | Front lens radius in a 4f system fibre collection system |
| f2Aperture | float | - | 1.0 | Back lens radius in a 4f system fibre collection system |
| frontOffset | float | - | 0.0 | distance between position and the front lens in a 4f system fibre collection system |
| backOffset | float | - | focalLenght2 | distance between fibre and back lens in a 4f system fibre collection system |
| frontToPinSep | float | - | focalLength1 | distance between front lens and a pinhole aperture in a 4f system fibre collection system |
| pinToBackSep | float | - | focalLength2 | distance between back lens and a pinhole aperture in a 4f system fibre collection system |
| pinAperture | float | - | max(f1Aperture, f2Aperture) | pinhole radius/size in a 4f system fibre collection system |
| acceptAngle | float | - | 90.0 | acceptance angle above the optical axis in degrees of the fibre in a 4f system fibre collection system |
| coreDiameter | float | - | 0.01 | diameter of the fibre core of the fibre in a 4f collection system |
| p1 | float array size 3 | - | [-1.0, -1.0, -1.0] | Used by camera detector only to set location and size of source|
| p2 | float array size 3 | - | [2.0, 0.0, 0.0] | See above |
| p3 | float array size 3 | - | [0.0, 2.0, 0.0] | See above |
| layer | integer | - | 1 | layer to match SDF layer label |
| nbins | integer | - | 100 | Number of bins in detector |
| maxval | float | - | 100.0 | Maximum value to bin |
| historyFileName | string | - | "photPos.obj" | Name of output file of detected photons histories |
| trackHistory | boolean | - | false | If true record detected photons histories. !!!Does not work with openMP!!! |

## Output

| Parameter | Type | Default | Notes |
|:---------:|:----:|:-------:|:-----:|
| fluence | string | fluence.nrrd | Filename for fluence output |
| absorb | string | absorb.nrrd | Filename for energy absorbed output |
| render_geomerty_name | string | geom_render.nrrd | Filename for render geometry output |
| render_geomerty | boolean | false | Render geometry out. For debugging purposes |
| render_source_name | string | source_render.nrrd | Filename for render source emission locations output |
| render_source | boolean | false | Render source emision locations out. For debugging purposes |
| render_size | integer array size 3 | [200, 200, 200] | Size in voxels of render |
| overwrite | boolean | false | Overwrite files if they have the same name |

## Simulation

| Parameter | Type | Default | Notes |
|:---------:|:----:|:-------:|:-----:|
| iseed | integer | 123456789 | seed for simulation. Each thread get its own copy + threadID |
| tev | boolean | false | Enables TEV image viewer to display simulation as it runs. Must have opened TEV prior to launching simulation. |
| absorb | boolean | false | Enables writing to file of absorbed energy. |

## Symmetry

| Parameter | Type | Options | Default | Notes |
|:---------:|:----:|:-------:|:-------:|:-----:|
| symmetryType | string | none, prism, flipped, uniformSlab, noneRotational, 360rotational | none | - |
| GridSize | float array size 3 | - | [10, 10, 10] | Number of voxels in symmetry grid in [x, y, z] for cartesian symmetry types, and [radius, theta, z] for rotational symmetry types |
! maxValues | float array size 3 | - | [1.0, 1.0 or 360 degrees, 1.0] | - |
| position | float array size 3 | - | [0.0, 0.0, 0.0] | Origin of symmetry grid with respect to absorption grid |
| direction | float array size 3 | - | [0.0, 0.0, 1.0] | Direction of z axis in symmetry grid with respect to absorption grid |
| rotation | float | 0.0 <= rotation < 360.0 | 0.0 | Rotation of x and y axis around the z axis |