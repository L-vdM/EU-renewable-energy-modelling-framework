# EU-renewable-energy-modelling-framework
This is a renewable energy modelling framework that was used to define meterological drivers of extreme events in the European renewable electrcity system, although it might be usefull for other puproses. It contains modules for the computatation of 1) wind energy production, 2) photvoltaic solar energy production, 3) run-of-river hydropower production, 4) reservoir hydropower inflow, 5) reservoir hydropower dispatch and 6) national demand. So far it has been run on daily resolution, but the model could be adjusted to deal with different timesteps. The current project is set-up to run at ERA5 grid for PV solar and wind production, and at 0.5x0.5 grid for run-of-river and reservoir hydropower prodcution.

This model has been validated against [ENTSO-E transparancy](https://transparency.entsoe.eu/) electricity production and demand data, and has been submitted for publication (April 2022). 


## Getting Started With This Framework

To get running, download this repository using the download button.
Create a new virtual environment and install the project requirements from the requirements.txt file.


### Define input data 

Then, set-up a new project for your own by:
- downloading (daily) meteorological input data (e.g. [ERA5](https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-pressure-levels?tab=overview)). Make sure to have a seperate netcdf for each variable. 
- regridding runoff data to the 0.5x0.5 grid of routing scheme (input_files/local-drain-direction/lddsound.nc)
- updating the references, directories and file names of the input data in the config.py file
- costumizing the input parameters as wanted in the config.py file

The framework can take seperate runs (e.g. a seperate netcdf file for each year, or an ensemble). It will use the climate files with the name of the run in it for each defined run. So make sure that netcdffiles that are input for the same run have the name of the run in the filename (could be a year).
NOTE: The hydropower dispatch model assumes perfect foresight for the length of the timeseries in each run. To take a longer period, timemerge the input files into a single run. 

If meteorological input data is is at a different grid than the installed capacities, make sure to  regrid the one of the two. 

The required climate input data is:

| Value                   | Description |
| :---                    | --- |
| tas                     | Daily mean near surface temperature [K]| 
| tasmax                  | Daily max near surface temperature [K]| 
| wind10m                 | Daily mean wind speed at 10 m above surface [m/s]| 
| runoff                  | Runoff  daily mean [m]|
| ssrd                    | Surface solar radiation downwards [W/m2]|

If you are uncertain about what to enter for any value or want to recreate the study as submitted for publication, just use the defaults.

You are now ready to get started,.

### Running the model
To run the model run the 00_run.py file.

## Generated Project Contents
Depending upon configurations when creating the project, the structure will look similar to the below:
The results can be found in the  output file. 


```                                                                                     
├── .00_run.py                                    <-  run this file to start run                  
├── agg_countries.py                              <-  module to aggregate countries, imported into 00_run.py.                     
├── compute_hydro_dispatch.py                     <-  module to acompute hydro dispatch, imported into 00_run.py.                                             
├── config.py                                     <-  configuration of the model run.
├── input_files                                   <-  direcotry with input files for run                  
│   ├── EEZ_land_union_v3_202003        <-  country borders                                             
│   ├── climate_data                    <-  directory for climate input data                                    
│   ├── demand_fit                      <-  directory with demand fit files                             
│   │   ├── demand_fit_values.nc                    <- the demand fit values per country        
│   │   ├── dict_population_per_country.json        <- dictory of country names and numbers
│   │   ├── population_t2m_grid.nc                  <- population per gridcell
│   │   └── population_t2m_grid_weights.nc          <- population weights of gridcell                                                   
│   ├── hydropower_dispatch             <- folder with hydropower dispatch input files                                          
│   │   ├── reservoir_filling_level_at_start.csv    <- filling level of reservoirs at start
│   │   └── storage_capacities.csv                  <- national storage capacities hydropower reservoirs (in weeks)
│   ├── installed-capacities            <- installed capacities 
│   │   ├── PVroof_cap_dunnett.nc                   <- PV installed capacity (default =0)                                                   
│   │   ├── PVutil_cap_dunnett.nc                   <- PV installed capacity                                    
│   │   ├── RES_reservoirHydro_cap_distr.nc         <- reservoir hydropower installed capacities                                        
│   │   ├── RES_run-of-riverHydro_cap_distr.nc      <- run-of-reservoir installed capacities                                    
│   │   ├── windOffhore_cap_MODnet.nc               <- offshore wind installed capacities                               
│   │   └── windOnshore_cap_dunnett.nc              <- onshore wind installed capacities                                    
│   ├── local-drain-direction                                                           
│   │   └── lddsound.nc       <- local drain direction file for routing runoff                                
│   └── mean_cf.csv                     <- mean national capacity factor of hydropower reservoirs                       
├── model                                         <- scripts used in model        
│   ├── constants                       <- mappings and constants used in modelling framework                           
│   │   ├── constant.py                                                             
│   │   ├── ldd.py                                                              
│   │   └── mappings.py                                                                 
│   ├── data_processing                 <-  scripts used for dataprocessing in modelling framework                                          
│   │   ├── attributes.py                                                               
│   │   ├── general.py                                                              
│   │   ├── masking.py                                                                  
│   │   └── production.py                                                               
│   ├── energy_computation              <-  scripts used in energy modules of framework                                             
│   │   ├── demand.py                                                       
│   │   ├── hydro.py                                                            
│   │   ├── solarpv.py                                                          
│   │   └── wind.py                                                         
│   ├── optimize_reservoirs             <- scripts used to generate reservoir dispatch                                      
│   │   ├── read_data.py                                                        
│   │   └── reservoir_model.py                                                      
│   └── routing                         <- script for routing                                                           
│       └── routing.py                                                 
├── ouput                                                          
│   └── run                             <- ouput directory. Set in config.                                  
├── readme.md                                     <- The README for people that would like the use this framework.                
├── requirements.txt                              <- The requirements file for reproducing the environment                    
└── route_runoff.py                               <- module to route runoff, imported into 00_run.py.
```

## Example input and output data
If you are interested in the output data from the work submitted to publish or have difficulty finding input data feel free to contact me.


## Other questions or Contributions to this framework
Contributions to this framework are appreciated and encouraged.

To contribute an update or if you have any other questions feel free to contact me.
