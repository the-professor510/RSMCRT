# Config file settings

The below sections describe the tables (dictionaries) that are able to be defined for SignedMCRT

## source

This table defines the parameters for the light source used in the simulation it can have the following:

| Parameter | Type | Options | Default | Notes |
|:---------:|:----:|:-------:|:-------:|:----:|
| name | string | point, circular, uniform, pencil, annulus, focus | point | - |
| nphotons | integer | - | 1000000 | - |
| position | float array size 3 | - | [0.0, 0.0, 0.0] | Default value only set for point source type|
| direction | float array size 3 or string | - | -z | String type applies to all source types bar: Uniform and circular |
| point1 | float array size 3 | - | [-1.0, -1.0, -1.0] | Used by uniform source only to set location and size of source |
| point2 | float array size 3 | - | [2.0, 0.0, 0.0] | See Above |
| point3 | float array size 3 | - | [0.0, 2.0, 0.0] | See Above |
| Radius | float | - | 0.5 | Only used by circular source |

## grid

| Parameter | Type | Default | Notes |
|:---------:|:----:|:-------:|:-----:|
| nxg | integer | 200 | Number of voxel in x direction |
| nyg | integer | 200 | Number of voxel in y direction |
| nyg | integer | 200 | Number of voxel in z direction |
| xmax | float | 1.0 | Half size of simulated medium in x direction |
| ymax | float | 1.0 | Half size of simulated medium in y direction |
| zmax | float | 1.0 | Half size of simulated medium in z direction |
| units | string | cm | Units of simulation (currently need to manually adjust optical properties to account) |

## geometry

| Parameter | Type | Default | Notes |
|:---------:|:----:|:-------:|:-----:|
| geom_name | string | sphere | Name of experiment for metadata |
| tau | float | 10.0 | Tau value for MCRT scattering test experiment |
| num_spheres | integer | 10 | Number of random spheres for sphere scene |
| musb | float | 0.0 | Optical properties for experimental geometry for whiskey Raman sensing paper |
| muab | float | 0.01 | See above |
| musc | float | 0.0 | See Above |
| muac | float | 0.01 | See Above |
| hgg | float | 0.7 | See Above |


## detectors

| Parameter | Type | Options | Default | Notes |
|:---------:|:----:|:-------:|:-------:|:----:|
| type | string | annulus, circle | - | - |
| position | float array size 3 | - | - | Central position of detector |
| radius1 | float | - | - | Radius of circular detector. Inner radius of annular detector |
| radius2 | float | - | - | Outer radius of annulus detector. Must be larger than radius1 |
| layer | integer | - | 1 | layer to match SDF layer label |
| nbins | integer | - | 100 | Number of bins in detector |
| maxval | float | - | 100.0 | Maximum value to bin |

## output

| Parameter | Type | Default | Notes |
|:---------:|:----:|:-------:|:-----:|
| fluence | string | fluence.nrrd | Filename for fluence output |
| render | string | geom_render.nrrd | Filename for render geometry output |
| render_geom | boolean | false | Render geometry out. For debugging purposes |
| render_size | integer array size 3 | [200, 200, 200] | Size in voxels of render |
| overwrite | boolean | false | Overwrite files if they have the same name |

## simulation

| Parameter | Type | Default | Notes |
|:---------:|:----:|:-------:|:-----:|
| iseed | integer | 123456789 | seed for simulation. Each thread get its own copy + threadID |
| tev | boolean | false | Enables TEV image viewer to display simulation as it runs. Must have opened TEV prior to launching simulation. |