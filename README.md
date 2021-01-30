# The Water Table Model (WTM)

***This model combines groundwater and surface water flow to output the elevation of the water table relative to the land surface at a given time.***

The model is intended for determining the depth or elevation of the water table, given a certain topography and set of climate inputs. Water table can be below ground (groundwater) or above ground (lake surfaces).

The model works by coupling groundwater and surface water components. The groundwater component moves water cell-to-cell using Darcy's Law with a finite-difference approach in a single layer of vertically integrated hydraulic conductivity. 

The surface-water component was collaboratively written by R Barnes and KL Callaghan. It works by creating a hierarchy of depressions for the topography, and then allowing water to move across the land surface, filling depressions and spilling from one depression into another. For more details on the depression hierarchy, see:

**Barnes, R, Callaghan, KL, and Wickert, AD, (2020), [Computing water flow through complex landscapes, Part 2: Finding hierarchies in depressions and morphological segmentations](https://esurf.copernicus.org/articles/8/431/2020/), *Earth Surf. Dynam.*, doi:10.5194/esurf-8-431-2020**

More details on the surface-water component, Fill-Spill-Merge, are available at:
**Barnes, R, Callaghan, KL, and Wickert, AD, (2020), [Computing water flow through complex landscapes, Part 3: Fill-Spill-Merge: Flow routing in depression hierarchies](https://esurf.copernicus.org/preprints/esurf-2020-31/), *Earth Surf. Dynam. Discuss.*, doi:10.5194/esurf-2020-31 **

This code has not been tested on Windows and may only work on Unix-based systems.

Please contact us if you have questions or suggestions!

## Required data inputs

Data inputs are all in a .nc (NetCDF) format. The NetCDF files should have 3 variables: 'lat' for latitude, 'lon' for longitude, and 'value' for the value of interest (e.g. elevation, precipitation, etc).

The following files are required:
* Topography - elevation in metres
* Slope - determined from your topography
* Mask - indicating the location of land (1) and ocean (0)
* Precipitation - in metres per year
* Evapotranspiration - in metres per year
* Winter air temperature - in degrees Celsius
* Hydraulic conductivity - in metres per second
* Porosity - unitless
* Open-water evaporation - in metres per second

## Dependencies

* The C++ compiler g++
* NetCDF for C++
* RichDEM

## Downloading with dependencies

The best way to obtain this code is by cloning the repository.

Before starting, note that in order to include RichDEM from GitHub, you will need a Public Key associated with your account. Instructions to do so can be found here.
https://help.github.com/en/github/authenticating-to-github/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent

Clone with submodule dependencies (RichDEM) included:
```sh
git clone --recurse-submodules https://github.com/KCallaghan/Coupled_SW_GW
```

If you forget to do this and just run a normal `git clone`, you can still pull the submodules:
```
git submodule update --init --recursive
```

In either case, use the following to update the submodules:
```
git pull --recurse-submodules
```

## Compilation
To build with `cmake` use:
```
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo ..
make
```
Use `-DSANITIZE_ADDRESS=On` to enable addressing sanitizing

## Running the code
Ensure that all of the data files are located appropriately in a folder together. Edit the Config_file.cfg configuration file as appropriate. The configuration file contains the following variables:

* textfilename       {The name of your output text file.txt}
* outfilename        {The name of your output depth-to-water-table file in geoTiff format.tif}
* cells_per_degree   {How many cells per degree in your data. E.g. One-degree resolution will be 1. 5 arcsecond resolution will be 12.}

time_start, surfdatadir, and region are all to help you name your input files or place them in a specific folder. The code will look for data files in the following format:
surfdatadir + region + time_start + "\_suffix.nc",
where the suffix refers to the specific file (topography, mask, precipitation, evaporation, winter_temperature, slope, open_water_evaporation, porosity, or ksat).
An example of a file path would be: "surfdata/North_America_10000_topography.nc".
In this case, you would set:

* surfdatadir        surfdata/
* region             North_America_
* time_start         10000

Two run types are possible: equilibrium, and transient. An equilibrium run assumes that the topography and climate are not changing and runs for many iterations until the equilibrium condition for the water table is found. Set the run_type parameter to 'equilibrium'.
A transient run requires a starting depth to water table as an additional input. The algorithm will then run for a set number of iterations, to represent a number of years passing, and output the new water table under a changing set of climatic and topographic conditions. In this case, both start and end states are required for all file inputs. Set the time_end parameter to lead to the files at the end time of the transient run, while time_start leads to the files at the initial time of the transient run. The run_type parameter should be set to 'transient'.
Other parameters include:

* deltat             {Number of seconds per time step, e.g. 315360000 for a 10-year time step}
* southern_edge      {Southern-most latitude of your domain in decimal degrees, e.g. 5}
* maxiter            500                  {How many times GW should run before FSM runs}
* total_cycles       1000                 {how many times FSM should run before completion} 
* infiltration_on    0                    {true is 1, false is 0. Only recommend true for high-resolution input data.}
* fdepth_a           200                  {e-folding depth coefficients}
* fdepth_b           150                  {e-folding depth coefficients}
* fdepth_fmin        2                    {e-folding depth coefficients}

Once the configuration file has been set up appropriately, simply open a terminal and type
```
./twsm.x Config_file.cfg
```
There will be some on-screen outputs to indicate the first steps through the code, after which values of interest will be output to the text file and an updated geoTiff output file will be saved every 100 iterations.

## Outputs
The program outputs a text file that provides information on the current minimum and maximum water table elevation, the changes in surface water and groundwater within the past iteration, and the number of iterations passed.
The main output is a geoTiff file that supplies the depth to/elevation of the water table. Negative values indicate a water table below the surface, while positive values indicate a water table above the surface (i.e. a lake).

## Completing a model run
A satisfactory method of detecting whether the model has reached equilibrium is still under construction. For now, it is at the discretion of the user whether he output after a given number of iterations is appropriate to use. The code will automatically complete after the number of iterations selected in the total_cycles parameter have been performed.
