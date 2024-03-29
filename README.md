[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4265369.svg)](https://doi.org/10.5281/zenodo.4265369)

# The Water Table Model (WTM)

***This model combines groundwater flow and dynamic lake simulation to output the elevation of the water table relative to the land surface at a given time.***

The model is intended for determining the depth or elevation of the water table, given a certain topography and set of climate inputs. Water table can be below ground (groundwater) or above ground (lake surfaces).

The model works by coupling groundwater and dynamic lake components. The groundwater component moves water cell-to-cell by solving the 2D groundwater flow equation for a heterogeneous, horizontally isotropic medium. It uses PETSc's SNES component to do so. The model is structured with a single layer of vertically integrated hydraulic conductivity.

The dynamic lake component was collaboratively written by R Barnes and KL Callaghan. It works by creating a hierarchy of depressions for the topography, and then allowing water to move across the land surface, filling depressions and spilling from one depression into another. For more details on the depression hierarchy, see:

**Barnes, R, Callaghan, KL, and Wickert, AD, (2020), [Computing water flow through complex landscapes, Part 2: Finding hierarchies in depressions and morphological segmentations](https://esurf.copernicus.org/articles/8/431/2020/), *Earth Surf. Dynam.*, doi:10.5194/esurf-8-431-2020**

More details on the surface-water component, Fill-Spill-Merge, are available at:
**Barnes, R, Callaghan, KL, and Wickert, AD, (2020), [Computing water flow through complex landscapes, Part 3: Fill-Spill-Merge: Flow routing in depression hierarchies](https://esurf.copernicus.org/preprints/esurf-2020-31/), *Earth Surf. Dynam. Discuss.*, doi:10.5194/esurf-2020-31 **

This code has not been tested on Windows and may only work on Unix-based systems.

Please contact us if you have questions or suggestions!

## Required data inputs

Data inputs are all in a Geotiff (.tif) format. 

The following files are required:
* Topography - elevation in metres
* Slope - determined from your topography
* Mask - indicating the location of land (1) and ocean (0)
* Precipitation - in metres per year
* Evapotranspiration - in metres per year
* Winter air temperature - in degrees Celsius
* Hydraulic conductivity - in metres per second
* Porosity - unitless
* Open-water evaporation - in metres per year

## Dependencies

* The C++ compiler g++
* GDAL
* RichDEM

## Downloading with dependencies

The best way to obtain this code is by cloning the repository.

Before starting, note that in order to include RichDEM from GitHub, you will need a Public Key associated with your account. Instructions to do so can be found here.
https://help.github.com/en/github/authenticating-to-github/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent

Clone with submodule dependencies (RichDEM) included:
```sh
git clone --recurse-submodules https://github.com/KCallaghan/WTM
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
cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo -DUSE_GDAL=ON ..
make
```
Use `-DSANITIZE_ADDRESS=On` to enable addressing sanitizing.

Alternatively, to build with `ninja`, use:
```
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo -DRICHDEM_LOGGING=ON  -GNinja ..
ninja
```

## Running the code
Ensure that all of the data files are located appropriately in a folder together. Edit the Config_file.cfg configuration file as appropriate. The configuration file contains the following variables:

* textfilename       {The name of your output text file.txt}
* outfile_prefix     {The name of your output depth-to-water-table file in geoTiff format. The code will append the time passed and .tif extension.}
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
# Optionally:
# export OMP_NUM_THREADS=N
# Required
./build/wtm.x Config_file.cfg -snes_mf -snes_type anderson -snes_stol T
```
Here, N is the number of CPU threads you want the parallel processing for the groundwater-flow step to use. In the above line, you are setting an environment variable that will define this until you exit the terminal window.
T is the tolerance setting for the SNES solver (e.g. 0.000001).

There will be some on-screen outputs to indicate the first steps through the code, after which values of interest will be output to the text file and an updated geoTiff output file will be saved every X iterations (X is set in the configuration file).

## Example of a full config file:
```
textfilename       my_model_run.txt     #The name of the output textfile, which will include printed values describing change in the water table.
outfile_prefix     my_model_run_        #An output file after 10 time steps will be named "my_model_run_0010.tif".
cells_per_degree   60                   #how many cells in one degree. This example for 1 arcminute cells.
#the below parameters are used in import file names. The code searches for files in the format:
#surfdatadir + region + time_start + input_type.tif,
#where the input_type is each of the input files discussed above (e.g. topography, precipitation, etc).
#For equilibrium runs, only time_start is used. For transient runs, time_start and time_end indicate
#the two sets of files with data at the beginning and at the end of the time period for the transient run.
time_start         021000               #used in filenames for import
time_end           020000               #used in filenames for import
surfdatadir        surfdata/            #used in filenames for import
region             my_area_             #used in filenames for import
deltat             31536000             #seconds in your timestep. This example for 1 year.
run_type           equilibrium          #test, equilibrium or transient
southern_edge      -52                  #Southern-most latitude of your import files
maxiter            1                    #how many times GW should run before FSM runs. Optionally, the groundwater can move multiple times before running FSM for
surface water. This is to save time on computation.
total_cycles       5000                 #how many times FSM should run before completion.
fdepth_a           100                  #e-folding depth coefficients
fdepth_b           150                  #e-folding depth coefficients
fdepth_fmin        2                    #e-folding depth coefficients
#should water be allowed to infiltrate during overland flow in FSM?
infiltration_on    0                    #true is 1, false is 0. Only recommend true for high-resolution input data.
#Are you supplying a starting water table? Note that you MUST supply a starting water table for transient runs.
#For equilibrium runs, you can optionally supply a starting water table; if you do not, water table will initialise = 0 everywhere.
supplied_wt        1                    #1 if you are supplying a starting water table, 0 if not
#Should surface water be moved using FSM? Only recommend turning this off for testing purposes.
fsm_on             1                    # 1 to enable Fill-Spill-Merge for routing surface water is enabled; 0 otherwise.
#Is water allowed to gather in lakes, with lake evaporation removing some portion of it?
#If this is set to 0, all surface water will be removed from the domain.
evap_mode          1                    # 1 to use a grid of potential evaporation for lakes; 0 to remove all surface water.
```

## Outputs
The program outputs a text file that provides information on the current minimum and maximum water table elevation, the changes in surface water and groundwater within the past iteration, and the number of iterations passed.
The main output is a geoTiff file that supplies the depth to/elevation of the water table. Negative values indicate a water table below the surface, while positive values indicate a water table above the surface (i.e. a lake).

## Completing a model run
A satisfactory method of detecting whether the model has reached equilibrium is still under construction. For now, it is at the discretion of the user whether he output after a given number of iterations is appropriate to use. The code will automatically complete after the number of iterations selected in the total_cycles parameter have been performed.
