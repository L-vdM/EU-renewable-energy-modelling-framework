import numpy as np
import glob

# =============================================================================
# 0. general sittings
# =============================================================================
dir_main = '/Users/lieke/Documents/01-diy-model/cleaned_for_online/'
# set directory of model function files
dir_functions = dir_main + 'model/'
# timestep
dt = 24
# set namee of author and project for output files 
author_name = 'author name'
project_name = 'project name'


# =============================================================================
# 1. set directories
# =============================================================================
## workstation
## computer
dir_rawdata = dir_main + 'input_files/'
dir_climate = dir_main + 'input_files/climate_data/'
# set output dir to create
dir_out = dir_main +'ouput/run3/'


# =============================================================================
# 2. set regions
# =============================================================================
# Set the countries; USE ISO31661A3 IDs to select countries
# or set 'europe'
selected_countries = 'europe'
## set_area in lon, lat 
l = -31.25
r = 74.5
t = 79
b = 33

# =============================================================================
# 3. input files
# =============================================================================
### initial boundary conditions
# general
# shapefile_countries = dir_rawdata + 'countries-shape/ne_10m_admin_0_countries.shp'
shapefile_countries = dir_rawdata + 'EEZ_land_union_v3_202003/EEZ_Land_v3_202030.shp'
country_name_column ='TERRITORY1'
abbreviation_column = 'ISO_TER1'

#for demand
population = dir_rawdata +'demand_fit/population_t2m_grid.nc'
pop_per_country = dir_rawdata +'demand_fit/population_per_country.nc'
demand_fit = dir_rawdata +'demand_fit/demand_fit_values.nc'

# routing
ldd_file = dir_rawdata + 'local-drain-direction/lddsound.nc'


# installed_capacities
# hydro_cap_csv = dir_rawdata + 'JRC-PPDB-OPEN.ver1.0/hydro_jrc.csv'
ror_cap = dir_rawdata + 'installed-capacities/RES_run-of-riverHydro_cap_distr.nc'
hydro_cap = dir_rawdata + 'installed-capacities/RES_ReservoirHydro_cap_distr.nc'
pv_util_cap = dir_rawdata + 'installed-capacities/PVutil_cap_dunnett.nc'
pv_roof_cap = dir_rawdata + 'installed-capacities/PVroof_cap_dunnett.nc'
onwind_cap = dir_rawdata + 'installed-capacities/windOnshore_cap_dunnett.nc'
offwind_cap = dir_rawdata + 'installed-capacities/windOffhore_cap_MODnet.nc'

### climate
runs = ['2020']
# 2m surface temperature
t2m_dir = dir_climate + 't2m_d/'#'tas_d/'
t2m_files = glob.glob(t2m_dir+"*.nc")
t2m_varname = 't2m'
# runoff
ro_dir = dir_climate + 'ro_d/'
ro_files = glob.glob(ro_dir+"*.nc")
ro_varname = 'ro'
# surface wind
wind_dir = dir_climate + 'wind10m_d/'#'sfcwind_d/'
wind_files = glob.glob(wind_dir+"*.nc")
wind_varname = 'wind10m'
# max daily temperature
t2mmax_dir = dir_climate + 't2mmax_d/'#'tasmax_d/'
t2mmax_files = glob.glob(t2mmax_dir+"*.nc")
t2mmax_varname = 't2mmax'
# radiation 
rad_dir = dir_climate + 'rsds_d/'
rad_files = glob.glob(rad_dir+"*.nc")
rad_varname = 'rsds'

# discharge
dis_dir = dir_out + 'discharge/'
dis_files = glob.glob(dis_dir+"*.nc")
dis_varname = 'discharge'



# =============================================================================
# 4. settings for routing - takes > 12 hours for 2000 years! 
# =============================================================================
route_runoff = True
route_runoff_submodules = dict(make_lddRoutingLayers = True)

# =============================================================================
# 5. settings for hydropower
# =============================================================================

# mean discharges to determine fraction of discharge used
dis_dir0 = dir_out + 'discharge/'
make_mean_yearly_discharge = True
mean_discharges_dir = dis_dir0
mean_discharge_file = dir_out + 'proccessed_bounds/mean_yearly_discharge.nc'

reservoir_efficiency = 0.9
ror_efficiency = 0.9

annual_hydro_cf = dir_rawdata + 'mean_cf.csv'

## q = 1- exceeding probability
q_ror = 0.75 

## reccession coefficient
kx = 0.6
# =============================================================================
# 5. settings for solar
# =============================================================================
# solar PV
cT_c1 = 4.3     # constant [dC]
cT_c2 = 0.943   # constant [-]
cT_c3 = 0.028   # constant [dC m2 W-1]
cT_c4 = -1.528  # constant [dC s m-1]
gamma = -0.005  # constant [--]
temp_ref = 25      # reference temperature [dC]

gstc = 1000     # standard test conditions [W m-2]
shift_doy = 186 # if HadGEM : 180
pv_constants = [cT_c1, cT_c2, cT_c3, cT_c4]

# =============================================================================
# 5. settings for wind
# =============================================================================
# wind
height_ref = 10.0   # height of available wind data [m]
time_oper = 24          # operational time of hub [h/day]


v_cutinland=2
v_ratedland=10
v_cutoutland=25

v_cutinoff = 3
v_ratedoff = 15
v_cutoutoff = 27

# variables
height_onshore = 90.0    # height onshore wind turbines [m]
height_offshore = 110.0 # height offshore wind turbines [m]
a_onshore = 0.143      # onshore roughness energy_type [-]
a_offshore = 0.11      # offshore roughness energy_type [-]

# =============================================================================
# 5. settings for optimization
# =============================================================================
storage_capacities = dir_rawdata + 'hydropower_dispatch/storage_capacities.csv' #country specific storage capacities
filling_level_at_start = dir_rawdata + 'hydropower_dispatch/reservoir_filling_level_at_start.csv' # initial 
default_stor_cap = '3M'## default cap storage of no cap storage in settings

# =============================================================================
# 6. settings for aggregate
# =============================================================================
agg = True

# =============================================================================
# 7. printing settings
# =============================================================================
#

def print_header(text):
        print(
    """
    -----------------------------------------------------
     %(h)s
    -----------------------------------------------------
    """ % {'h': text}
     )