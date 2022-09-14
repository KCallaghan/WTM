[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4265369.svg)](https://doi.org/10.5281/zenodo.4265369)

# The Water Table Model (WTM)

***This model combines groundwater and surface water flow to output the elevation of the water table relative to the land surface at a given time.***

The model is intended for determining the depth or elevation of the water table, given a certain topography and set of climate inputs. Water table can be below ground (groundwater) or above ground (lake surfaces).

The model works by coupling groundwater and surface water components. The groundwater component moves water by solving the 2D horizontal groundwater equation: 

<img src="https://latex.codecogs.com/svg.image?S_0\frac{\partial&space;h}{\partial&space;t}&space;=&space;\frac{\partial}{\partial&space;x}&space;\left(T\frac{\partial&space;h}{\partial&space;x}\right)&space;&plus;&space;&space;\frac{\partial}{\partial&space;y}&space;\left(T\frac{\partial&space;h}{\partial&space;y}\right)&space;&plus;&space;R" title="https://latex.codecogs.com/svg.image?S_0\frac{\partial h}{\partial t} = \frac{\partial}{\partial x} \left(T\frac{\partial h}{\partial x}\right) + \frac{\partial}{\partial y} \left(T\frac{\partial h}{\partial y}\right) + R" />

Here, h is the groundwater head, T is transmissivity (and is dependent on h), t is the length of a single user-defined time step, x and y are the two dimensions of groundwater movement, R is groundwater recharge, and S_0 is the specific storage through the thickness of the aquifer (which may be approximated as being equal to porosity). We use a single layer of vertically integrated hydraulic conductivity, and model all groundwater flow as saturated flow. 

The surface-water component was collaboratively written by R Barnes and KL Callaghan. It works by creating a hierarchy of depressions for the topography, and then allowing water to move across the land surface, filling depressions and spilling from one depression into another. For more details on the depression hierarchy, see:

**Barnes, R, Callaghan, KL, and Wickert, AD, (2020), [Computing water flow through complex landscapes, Part 2: Finding hierarchies in depressions and morphological segmentations](https://esurf.copernicus.org/articles/8/431/2020/), *Earth Surf. Dynam.*, doi:10.5194/esurf-8-431-2020**

More details on the surface-water component, Fill-Spill-Merge, are available at:
**Barnes, R, Callaghan, KL, and Wickert, AD, (2020), [Computing water flow through complex landscapes, Part 3: Fill-Spill-Merge: Flow routing in depression hierarchies](https://esurf.copernicus.org/preprints/esurf-2020-31/), *Earth Surf. Dynam. Discuss.*, doi:10.5194/esurf-2020-31 **

This code has not been tested on Windows and may only work on Unix-based systems.

Please contact us if you have questions or suggestions!

## Required data inputs

Data inputs are all in a geotiff (.tif) format. The following files are required:

* Topography - elevation in metres
* Slope - this should be determined from your topography
* Mask - indicating the location of land (1) and ocean (0). 
* Precipitation - in metres per year
* Evapotranspiration - in metres per year
* Winter air temperature - in degrees Celsius
* Horizontal hydraulic conductivity - in metres per second
* Porosity - unitless
* Open-water evaporation - in metres per year

Additionally, there are the following optional input files:

* Runoff ratio - unitless
* Starting water table - in metres relative to the land surface (i.e. topography). A starting water table is a required input for transient runs. 
* Vertical hydraulic conductivity - in metres per second 

Note that for transient model runs (discussed more below), the code will require separate input files for the start and the end times for topography, slope, precipitation, evapotranspiration, open water evaporation, winter temperature, and runoff ratio (if using).

## Dependencies

* Petsc [https://github.com/petsc/petsc]
* The C++ compiler g++
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
cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo -DUSE_GDAL=ON ..
make
```
Use `-DSANITIZE_ADDRESS=On` to enable addressing sanitizing

## Running the code
Ensure that all of the data files are located appropriately in a folder together. Edit the Config_file.cfg configuration file as appropriate. The configuration file contains the following variables:

* textfilename       {The name of your output text file.txt. This can be any name that you want, e.g. 'test.txt'. This file will contain text information about the progression of your WTM run, detailed more below.}
* outfile_prefix     {The name of your output depth-to-water-table file in geoTiff format. The code will append the time passed and .tif extension. For example, name this 'my_water_table_output_' and output files will take the format 'my_water_table_output_100.tif'.}
* cells_per_degree   {How many cells per degree in your data. E.g. One-degree resolution will be 1; 5 arcsecond resolution will be 12; 30 arcsecond resolution will be 120.}

* time_start
* time_end
* surfdatadir
* region

time_start, surfdatadir, and region are all to help you name your input files or place them in a specific folder. The code will look for data files in the following format:
surfdatadir + region + time_start + "\_suffix.tif",
where the suffix refers to the specific file (topography, mask, precipitation, evaporation, winter_temperature, slope, open_water_evaporation, porosity, horizontal_ksat, starting_wt, vertical_ksat, or runoff_ratio).
An example of a file path would be: "surfdata/North_America_10000_topography.tif".
In this case, you would set:

* surfdatadir        surfdata/
* region             North_America
* time_start         10000

* run_type

Three run types are possible: equilibrium, transient, and test.

**Test runs:**
Test runs are, as the name implies, only for testing. Use these as needed to help diagnose problems with your inputs or outputs. Test runs only require inputs for slope and topography, and use preset values for everything else. To use this, set the run_type parameter in the configuration file to 'test'.

**Equilibrium runs:**
An equilibrium run assumes that the topography and climate are not changing and complete many model time steps until the equilibrium (steady-state) condition for the water table is found. This scenario is most realistic for times and locations in Earth history where climate and topography has not changed significantly over the prior thousands of years. Because transient runs require a starting water table, it may also be used to determine this starting water table, but do note when doing so that you are assuming the water table at the start time was in equlibrium with climate and topography. To use this, set the run_type parameter in the configuration file to 'equilibrium'.

**Transient runs:**
A transient run requires a starting depth to water table as an additional input (this is optional for equilibrium runs). WTM will then run for a set number of time steps, and output the new water table under a changing set of climatic and topographic conditions. In this case, both start and end states are required for the relevant file inputs. Set the time_end parameter to lead to the files at the end time of the transient run, while time_start leads to the files at the initial time of the transient run. The run_type parameter should be set to 'transient'.

Continuing to the rest of the configuration file parameters:

* deltat             {Number of seconds per time step, e.g. 31536000 for a 1-year time step}
* southern_edge      {Southern-most latitude of your domain in decimal degrees, e.g. 5. This is used to compute cell sizes across the domain.}
* maxiter            10                   {How many times the groundwater portion of the code should run before FSM runs}
* total_cycles       1000                 {how many times FSM should run before before exiting the program. If deltat, above, were set to 1 year (i.e. 31536000 seconds), maxiter is set to 10, and total_cycles is set to 1000, then the total real-world time that the program will simulate is 1*10*1000 = 10000 years.}
* infiltration_on    0                    {true is 1, false is 0. Only recommend true for high-resolution input data.}
* fdepth_a           200                  {e-folding depth coefficients}
* fdepth_b           150                  {e-folding depth coefficients}
* fdepth_fmin        2                    {e-folding depth coefficients}
* supplied_wt        0                    {a flag to set whether or not a starting estimated water table is supplied, for equilibrium runs. 0 is false, i.e. no starting water table is supplied, and the program will use a default value. 1 is true, i.e. there is a starting water table file supplied.}
* runoff_ratio_on    0                    {a flag to set whether or not a runoff ratio file is supplied. 0 is false, i.e. no runoff ratio is supplied, and the program will assume all recharge moves directly to the water table. 1 is true, i.e. there is a runoff ratio file supplied.}
* fsm_on             1                    {a flag to set whether or not the program should run fill-spill-merge. In general, we recommend leaving this on, but this may depend on your needs. 1 means FSM will run; 0 means it will not run.}
* evap_mode          1                    {a flag to set whether or not the program should use the surface-water evaporation file you have supplied. If this is set to 0, then **all surface water will be removed from the system**. No lakes will be allowed to form, and this will ultimately affect the groundwater table as well. Set to 1 to use your input evaporation files; 0 to remove all surface water}
* cycles_to_save     100                  {set to any integer value. The program will save an output file each time it has completed this many time steps. It will also always save an output file after completing the full number of steps in total_cycles.}

Once the configuration file has been set up appropriately, simply open a terminal and type
```
# Optionally:
# export OMP_NUM_THREADS=N
# Required
./build/wtm.x Config_file.cfg -snes_mf -snes_type anderson -snes_stol M
```
Here, N is the number of CPU threads you want the parallel processing for the groundwater-flow step to use and M is the tolerance to use for petsc's SNES solver, for example, 0.000001.

There will be some on-screen outputs to indicate the first steps through the code, after which values of interest will be output to the text file and an updated geoTiff output file will be saved every cycles_to_save number of iterations.

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
region             my_area             #used in filenames for import
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
runoff_ratio_on    1                    # 1 to use the runoff ratio input; 0 to assume all water goes directly to the water table
cycles_to_save     100                  # frequency with which to save outputs
```

## Outputs
The program outputs a text file that provides information on the current minimum and maximum water table elevation, the changes in surface water and groundwater within the past iteration, and the number of iterations passed.
The main output is a geoTiff file that supplies the depth to/elevation of the water table. Negative values indicate a water table below the surface, while positive values indicate a water table above the surface (i.e. a lake).

## Completing a model run
A satisfactory method of detecting whether the model has reached equilibrium is still under construction. For now, it is at the discretion of the user whether he output after a given number of iterations is appropriate to use. The code will automatically complete after the number of iterations selected in the total_cycles parameter have been performed.
