# Coupled_SW_GW

***This model combines groundwater and surface water flow to output the elevation of the water table relative to the land surface at a given time.***

The model is intended for determining the depth or elevation of the water table, given a certain topography and set of climate inputs. Water table can be below ground (groundwater) or above ground (lake surfaces). 

The model works by coupling groundwater and surface water components. The groundwater component is based on the methods used by Fan et al in their paper:

**Fan, Y, Li, H, and Miguez-Macho, G, (2013),[Global Patterns of Groundwater Table Depth](https://science.sciencemag.org/content/339/6122/940.abstract),*Science*, doi:10.1126/science.1229881**

The surface-water component was collaboratively written by R Barnes and KL Callaghan. It works by creating a hierarchy of depressions for the topography, and then allowing water to move across the land surface, filling depressions and spilling from one depression into another. For more details on the depression hierarchy, see:

**Barnes, R, Callaghan, KL, and Wickert, AD, (2019), [Computing water flow through complex landscapes, Part 2: Finding hierarchies in depressions and morphological segmentations](https://www.earth-surf-dynam-discuss.net/esurf-2019-34/),*Earth Surf. Dynam. Discuss.*, doi:10.5194/esurf-2019-34**

More details on the working of surface-water movement across the land surface will be added soon. 

This code has not been tested on Windows and may only work on Unix-based systems. 

Please contact us if you have questions or suggestions! 

## Required data inputs

Data inputs are all in a .nc (NetCDF) format. The NetCDF files should have 3 variables: 'lat' for latitude, 'lon' for longitude, and 'value' for the value of interest (e.g. elevation, precipitation, etc). 

The following files are required:
Topography - elevation in metres
Mask - indicating the location of land (1) and ocean (0)
Precipitation - in metres per year
Evaporation - in metres per year
Relative humidity - as a proportion from 0 to 1
Winter temperature - in degrees Celsius
E-folding depth - see the Fan et al paper linked above for details on how this is computed using slope and temperature
Hydraulic conductivity - in metres per second

##Dependencies

* The C++ compiler g++
* NetCDF for C++
* RichDEM

## Compilation
Use the supplied Makefile. Open a terminal in the folder containing the makefile and source code, and type
```make```

## Running the code
Ensure that all of the data files are located appropriately in a folder together. Edit the global.cfg configuration file as appropriate. The configuration file contains the following variables:

textfilename       {The name of your output text file.txt}
outfilename        {The name of your output depth-to-water-table file in NetCDF format.nc}
cells_per_degree   {How many cells per degree in your data. E.g. One-degree resolution will be 1. 5 arcsecond resolution will be 12.}

time_start, surfdatadir, and region are all to help you name your input files or place them in a specific folder. The code will look for data files in the following format:
surfdatadir + region + time_start + "\_suffix.nc",
where the suffix refers to the specific file (topo, mask, precip, evap, relhum, temp, fdepth, or ksat). 
An example of a file path would be: "surfdata/North_America_10000_topo.nc".
In this case, you would set:

surfdatadir        surfdata/
region             North_America_
time_start         10000

time_end and run_type will become relevant in a future release, for now leave these as they are. 

deltat             {Number of seconds per time step, e.g. 315360000 for a 10-year time step}
southern_edge      {Southern-most latitude of your domain in decimal degrees, e.g. 5}
maxiter            {How many iterations groundwater should run before running surface water. If in doubt, leave this as 1 at first.}

Once the configuration file has been set up appropriately, simply open a terminal and type 
```
./a.out global.cfg
```
There will be some on-screen outputs to indicate the first steps through the code, after which values of interest will be output to the text file and an updated netcdf output file will be saved every 10 iterations. 

## Outputs
The program outputs a text file that provides information on the current minimum and maximum water table elevation, the changes in surface water and groundwater within the past iteration, and the number of iterations passed. 
The main output is a netcdf file that supplies the depth to/elevation of the water table. Negative values indicate a water table below the surface, while positive values indicate a water table above the surface (i.e. a lake). 

## Completing a model run
A satisfactory method of detecting whether the model has reached equilibrium is still under construction. For now, it is at the discretion of the user whether he output after a given number of iterations is appropriate to use. The code will automatically complete after 500000 iterations have been performed.
