## test demand production
import glob
import os
import sys
import xarray as xr
import numpy as np
import pandas as pd

from pathlib import Path
from cdo import Cdo
cdo = Cdo()

import config
import shutil

# set function path
sys.path.append(config.dir_functions)

# run pre-processing modules
if config.route_runoff is True:
    import route_runoff

import constants.mappings as mappings

# data processing
import data_processing.masking as masking
import data_processing.attributes as attributes

# energy computations
import energy_computation.demand as demand
import energy_computation.hydro as hydro2
import energy_computation.solarpv as solarpv
import energy_computation.wind as wind

# =============================================================================
# Set files
# =============================================================================
### import files
df_countries = mappings.df_countries_select
# demand uses EU_map so drop nan values
df_countries = df_countries.loc[df_countries.EU_map.notna()]
# boundary import files
global_population = config.population
shapefile_countries = config.shapefile_countries
# import pre_processed data
fitvalues_file = config.demand_fit
dis_files = glob.glob(config.dis_dir + "*.nc")

### set output files
# directories
dir_out = config.dir_out
dir_bounds = dir_out + "proccessed_bounds/"
dir_demand = dir_out + "energy_demand/"
dir_prod = dir_out + "energy_prod_wind_solar/"
dir_reservoir = dir_out + "energy_inflow_reservoirs/"
dir_ror = dir_out + "energy_prod_ror/"
# create output directories 
for d in [dir_bounds, dir_demand, dir_reservoir, dir_prod, dir_ror]:
    isExist = os.path.exists(d)
    if not isExist:
        # create a new directory when it does not exist
        os.makedirs(d)
# processed boundary files
population_tempgrid = config.population
# temporary files
population_per_country = config.pop_per_country
# country_name_index_mapping = dir_bounds + "/dict_" + os.path.basename(population_per_country)[:-3] +".json"

## save a copy of config file in runfolder
shutil.copyfile("config.py", dir_out+"config_of_run.py")
# =============================================================================
# Set vars
# =============================================================================
var1 = config.t2m_varname # temperature
var2 = config.dis_varname
dt = config.dt

# =============================================================================
# Compute runoffriver hydropower
# =============================================================================
config.print_header("COMPUTING RUNOFFRIVER HYDROPOWER")
# import capacities
ror_cap = config.ror_cap

### define factor of discharge in gridcell
# discharge in a mean year
dis_formean = xr.open_mfdataset(config.mean_discharges_dir + "*nc")
meanyear = dis_formean[config.dis_varname].groupby("time.dayofyear").mean().compute()
dis_formean.close()
# exceeding probability flow
q = config.q_ror
dis_ep = meanyear.quantile(q, dim="dayofyear")

# capacity and dicharge at max capacity
cap = xr.open_dataset(ror_cap)
dis_cap = hydro2._compute_m3_capacity(
    cap.capacity, cap.weighted_H, eta=config.ror_efficiency
	) * (3600 * dt)  # make daily values

# Assign factor of cell
f = xr.where((dis_cap / dis_ep) > 1, 1, (dis_cap / dis_ep))

for file in dis_files:
    dis = xr.open_dataset(file)
    production = hydro2.compute_hydro_energy_production(
        dis[config.dis_varname] * f, cap.weighted_H
    ).compute()
    prod = attributes.set_global_attributes(
            production.to_dataset(name="ror"), 
            "jrc database and ERA5", 
            grid="0.5degx05deg", 
            area="Europe")
    prod = prod.drop("quantile")
    prod.ror.attrs.update(
             standard_name = "ror production",
            long_name =  "daily hydropower run-of-river production",
            units = "MWh")
    cap = cap.sel(lat=prod.lat, lon=prod.lon)
    prod_capped = xr.where(prod>(cap.capacity * dt), (cap.capacity * dt), prod)
    
    prod_capped.to_netcdf(dir_ror + os.path.basename(file).replace("dis", f"ror{q}"))
    print("done with :" + os.path.basename(file).replace("dis", f"ror{q}"))


# =============================================================================
# Compute demand
# =============================================================================
config.print_header("COMPUTING DEMAND")

### make population weights
masking.cut_netcdf_into_regions(
    df_countries.EU_map,
    global_population,
    population_per_country,
    shapefile_countries,
    country_indexes=df_countries.index_nr,
)
# regrid to temperature grid
cdo.remapsum(
    config.t2m_files[0],
    input=population_per_country,
    output=population_tempgrid,
    readCdf=True,
    options="-f nc",
)
# take weights
pop_temp = xr.open_dataset(
    population_tempgrid,
)
weights = pop_temp / pop_temp.sum(dim=["lat", "lon"])
weights.to_netcdf(population_tempgrid[:-3] + "_weights.nc")
# remove temp files
os.remove(population_per_country)

### compute demand
fv = xr.open_dataset(fitvalues_file)
# loop over temperature files

for r in config.runs:
    print("Computing demand for", r)
    # get file names of run
    tempf = [f for f in config.t2m_files if r in f][0]
    # open temperature data
    ds_t2m = xr.open_dataset(tempf)
    # prepare temperature dataset
    try:
        if ds_t2m[var1].units == "K":
            ds_t2m[var1] = ds_t2m[var1] - 273.15
            ds_t2m[var1].attrs["units"] = "degC"
    # if no units assume Kelvin
    except:
        ds_t2m[var1] = ds_t2m[var1] - 273.15
        ds_t2m[var1].attrs["units"] = "degC"
    ds_t2m[var1].attrs.update(standard_name=var1)
    
    # split run over years to reduce memory usage
    demand_list = []
    for y in np.unique(ds_t2m[var1]["time.year"].values):
        print(y)
        ds_demand0 = xr.Dataset()
        # weighted temperature
        ds_demand0["temp"] = (ds_t2m[var1].sel(time=str(y)) * weights.population).sum(
            dim=["lat", "lon"], keep_attrs=True
        )
        demand_list.append(ds_demand0)
    # concat years together again
    ds_demand = xr.concat(demand_list, dim="time") 

    # to match dimensions of fitvalues (country,period) per weekend and weekday
    ds_demand = xr.concat(
        [
            ds_demand.where(ds_demand["time.dayofweek"] < 5, drop=True),
            ds_demand.where(ds_demand["time.dayofweek"] >= 5, drop=True),
        ],
        "period",
    )
    ds_demand["period"] = ["weekday", "weekend"]
    # select only the countries for which we have fitted data:
    # select only the countries for which we have fitted data and for which we have climate data:
    demand_countries = [f for f in fv.country.data if f in ds_demand.country.values]
    ds_demand = ds_demand.sel(country=demand_countries)
    fv = fv.sel(country=demand_countries)

    # compute demand with fit variables
    ds_demand = demand.compute_demand(ds_demand, fv)
    # remove period dimension
    ds_demand = xr.concat([ds_demand.isel(period=0), ds_demand.isel(period=1)], "time")
    # update attributes
    ds_demand = attributes.set_global_attributes(
        ds_demand, "Entsoe-ERA5 fit and HW3", grid="gaussian n80", area="Europe"
    )

    # clean file
    ds_demand["demand"] = ds_demand.demand.transpose("time", "country")
    ds_demand = ds_demand.drop("period").dropna(dim="time")
    # sort the weekend and weekday demand by time
    ds_demand = ds_demand.sortby("time")
    # save file
    ds_demand.to_netcdf(dir_demand + os.path.basename(tempf).replace(var1, "demand"))
    print("done with :" + os.path.basename(tempf).replace(var1, "demand"))


# =============================================================================
# Compute reservoir hydropower
# =============================================================================
config.print_header("COMPUTING RESERVOIR HYDROPOWER")

# import capacities
reservoir_cap = config.hydro_cap

# make mean yearly discharge file
if config.make_mean_yearly_discharge:
    cdo.timmean(
        input="-yearsum -ensmean %s" % config.mean_discharges_dir + "*.nc",
        output=config.mean_discharge_file,
    )
entsoe_eu = df_countries.entsoe_transparency.replace("", np.nan).dropna()
mappingcf = dict(zip(entsoe_eu.index, entsoe_eu))

### set variables
n = config.reservoir_efficiency
### open data
disy = xr.open_dataset(config.mean_discharge_file).isel(time=0)
res = xr.open_dataset(reservoir_cap)
res = masking.cut_box(res)
cf = pd.read_csv(config.annual_hydro_cf, index_col=0)

### add capacity factor to hydropower gridcells
# mean capacity factor over years available
cf = cf.mean()
cf.index = cf.index.map(mappingcf)
cf.name = "cf"
# add capacity factor to
gridcell_with_hydro = res.country_code.to_dataframe().reset_index()
gridcell_with_hydro.replace("", np.nan, inplace=True)
gridcell_with_hydro = gridcell_with_hydro.dropna().set_index("country_code")
cf_per_gridcell = gridcell_with_hydro.join(cf)
cf_per_gridcell = cf_per_gridcell.set_index(["lat", "lon"]).to_xarray()
res["cf_country"] = cf_per_gridcell.cf
res["dis_ratio"] = hydro2.determine_discharge_ratio(
    res.capacity, res.weighted_H, res.cf_country, disy[var2]
)

### compute production for every discharge file
# for r in config.runs:
#     print("Computing discharge for", r)
# # get file names of run
#     disf = [f for f in config.dis_files if r in f][0]

for disf in dis_files:
    dis = xr.open_dataset(disf)
    dis = masking.cut_box(dis)

    production = hydro2.compute_hydro_energy_production(
        dis[var2] * res["dis_ratio"], res.weighted_H
    )
    production = production.to_dataset(name="inflow")
    production.inflow.attrs.update(
        standard_name="Hydropower reservoir inflow",
        long_name="Hydropower reservoir inflow",
        units="MWh",
    )
    production.to_netcdf(
        dir_reservoir + os.path.basename(disf).replace(var2[:3], "reservoir_in")
    )
    print("done with :" + os.path.basename(disf).replace(var2[:3], "reservoir_in"))


# =============================================================================
# PV solar and wind production
# =============================================================================
config.print_header("COMPUTING SOLAR AND WIND PRODUCTION")

dir_bounds = config.dir_out + "proccessed_bounds/"
capacity_files = [
    config.pv_util_cap, 
    config.pv_roof_cap,
    config.onwind_cap,
    config.offwind_cap
]
# remap the capacities to temperature file
for cf in capacity_files:
    cdo.remapsum(
        config.t2m_files[0],
        input=cf,
        output=dir_bounds + Path(cf).stem + "_remap.nc",
        readCdf=True,
        options="-f nc",
    )
pv_util_cap = xr.open_dataset(dir_bounds + Path(config.pv_util_cap).stem + "_remap.nc")
pv_roof_cap = xr.open_dataset(dir_bounds + Path(config.pv_roof_cap).stem + "_remap.nc")
onwind_cap = xr.open_dataset(dir_bounds + Path(config.onwind_cap).stem + "_remap.nc")
offwind_cap = xr.open_dataset(dir_bounds + Path(config.offwind_cap).stem + "_remap.nc")


for r in config.runs:
    print("Computing PV and wind production for", r)
    prod = xr.Dataset()
    # get file names of run
    tempf = [f for f in config.t2m_files if r in f][0]
    tempmaxf = [f for f in config.t2mmax_files if r in f][0]
    sfcwindf = [f for f in config.wind_files if r in f][0]
    radf = [f for f in config.rad_files if r in f][0]
    # open datasets
    temp = xr.open_dataset(tempf)
    tempmax = xr.open_dataset(tempmaxf)
    sfcwind = xr.open_dataset(sfcwindf)
    rad = xr.open_dataset(radf)

    # solar
    pot_pv = solarpv.compute_solar_energy_potential(
        rad[config.rad_varname],
        temp[config.t2m_varname],
        tempmax[config.t2mmax_varname],
        sfcwind[config.wind_varname],
    )

    prod["pv_util"] = solarpv.compute_solar_energy_production(
        pot_pv, pv_util_cap.CAP, energy_type="utility solar"
    )

    prod["pv_roof"] = solarpv.compute_solar_energy_production(
        pot_pv, pv_roof_cap.CAP, energy_type="rooftop solar"
    )

    # wind
    pot_wind_off = wind.compute_wind_energy_potential(
        sfcwind[config.wind_varname],
        config.height_offshore,
        config.a_offshore,
        energy_type="offshore wind",
        v_cutin=config.v_cutinoff,
        v_rated=config.v_ratedoff,
        v_cutout=config.v_cutoutoff,    )

    pot_wind_on = wind.compute_wind_energy_potential(
        sfcwind[config.wind_varname],
        config.height_onshore,
        config.a_onshore,
        energy_type="onshore wind",
        v_cutin=config.v_cutinland,
        v_rated=config.v_ratedland,
        v_cutout=config.v_cutoutland,
    )

    prod["wind_offshore"] = wind.compute_wind_energy_production(
        pot_wind_off, offwind_cap.CAP, energy_type="offshore wind"
    )

    prod["wind_onshore"] = wind.compute_wind_energy_production(
        pot_wind_on, onwind_cap.CAP, energy_type="onshore wind"
    )

    prod = attributes.set_global_attributes(
        prod, f"climate data from run {r}, capacities_eu.nc"
    )

    prod.to_netcdf(
        dir_prod + os.path.basename(tempf).replace(config.t2m_varname, "prod")
    )

    print("done with :" + os.path.basename(tempf).replace(config.t2m_varname, "prod"))


if config.agg == True:
    import agg_countries
if config.optimization  == 'rolling':
    import compute_rolling_hydro_dispatch
elif config.optimization == 'foresight':
    import compute_hydro_dispatch
else:
    print('optimization type unknown. choose between rolling or foresight')