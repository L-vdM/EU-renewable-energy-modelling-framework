import pandas as pd
import pyomo.environ as pyo
# Create a solver
opt = pyo.SolverFactory('glpk')
import config

import numpy as np
import pandas as pd
# import matplotlib.pyplot as plt
import xarray as xr
import os 
import glob
import config
import sys
# set function path
sys.path.append(config.dir_functions)
import constants.mappings as mappings
import optimize_reservoirs.reservoir_model as rm
import pyomo.environ as pyo
from pyomo.opt import SolverFactory

# from scipy.stats import pearsonr
# from scipy.stats import linregress

import data_processing.attributes as attributes


## country mapping 
df_countries = mappings.df_countries_select


storage_capacityf = config.storage_capacities
fill_fracf = config.filling_level_at_start
best_storage = pd.read_csv(storage_capacityf, index_col=0)
best_storage.columns = ['week']
# best_storage = best_storage.replace(1, 4)
best_storage = best_storage.replace(np.nan, 1)
fill_frac = pd.read_csv(fill_fracf, index_col=0)

## Run optimization
# Create a solver
# opt = pyo.SolverFactory('ipopt')
opt = pyo.SolverFactory('glpk')


config.print_header(f'COMPUTING HYDRO OUTFLOW')

## open files
### import files
df_countries = mappings.df_countries_select
# demand uses EU_map so drop nan values 
df_countries = df_countries.loc[df_countries.EU_map.notna()]
# make mappings
mapping = dict(zip(df_countries.index_nr, df_countries.EU_map))
mapping_cap =dict(zip(df_countries.nuts_id, df_countries.index))
mapping_back = dict(zip(df_countries.EU_map, df_countries.index_nr))


# directories
dir_bounds = config.dir_out + 'proccessed_bounds/'
dir_agg = config.dir_out + 'agg_production/'
dir_demand = config.dir_out + 'energy_demand/'
dir_pv_wind = dir_agg + 'energy_prod_wind_solar/'
dir_res_in = dir_agg + 'energy_inflow_reservoirs/'
dir_ror = dir_agg + 'energy_prod_ror/'

### set output directories
dir_res_out =  dir_agg + 'energy_outflow_reservoirs/'
# make output directory
for d in [dir_res_out]:
    isExist = os.path.exists(d)
if not isExist:
    # create a new directory when it does not exist
    os.makedirs(d)


# open country capacities based on gridded input
reservoir_cap = config.hydro_cap
df_cap = xr.open_dataset(reservoir_cap).to_dataframe()
df_cap = df_cap.groupby('country_code').sum()
df_cap = df_cap[['capacity']]
# clean data 
df_cap.replace("", np.nan, inplace=True)
df_cap.replace(0, np.nan, inplace=True)
df_cap.dropna(inplace=True)
# map country_code
df_cap.index = [mapping_cap[k] for k in df_cap.index]

### open files
demand_files = glob.glob(dir_demand+"*.nc")
inflow_files = glob.glob(dir_res_in+"*.nc")
production_files = glob.glob(dir_pv_wind+"*.nc")
ror_prod_files = glob.glob(dir_ror+"*.nc")

for r in config.runs:
    print('Computing dispatch for', r)
    # set infiles 
    demand_file = [f for f in demand_files if r in f][0]
    inflow_file = [f for f in inflow_files if r in f][0]
    production_file = [f for f in production_files if r in f][0]
    ror_file =  [f for f in ror_prod_files if r in f][0]

    # demand
    demand = xr.open_dataset(demand_file)
    demand['country'] = list(map(mapping.get, demand['country'].values))
    demand = demand.dropna(dim='time')

    # inflow
    inflow = xr.open_dataset(inflow_file)
    inflow['country'] = list(map(mapping.get, inflow['country'].values))
    # inflow['time'] = [t + pd.Timedelta('+11hours30Min') for t in inflow['time'].values]
    # production
    prod_m = xr.open_dataset(production_file)
    prod_m['country'] = list(map(mapping.get, prod_m['country'].values))

    # run of river production
    ror_prod = xr.open_dataset(ror_file)
    ror_prod['country'] = list(map(mapping.get, ror_prod['country'].values))
    # sum all vars in production (wind and solar)
    ror_prod = ror_prod.loc[dict(time=ror_prod.time.where(demand.time))]
    prod_m['ror'] = ror_prod.ror

    # sum all vars in production (wind, solar and ror)
    vars_to_sum = list(prod_m.keys())
    prod_total = prod_m[vars_to_sum].to_array().sum("variable").to_dataset(name='wind_solar_ror')


    appended_data = []
    countries = []
    for c in demand.country.values:
    #     print(c)
        dc = demand.sel(country=c).sortby('time').to_dataframe()
        try:
            # open reservoir inflow and drop 0 values
            ic = inflow.sel(country=c).to_dataframe()
            ic = ic.replace(0, np.nan)
        except:
            print(f'{c} not in inflow file')
            # no need to optimize if country has no hydro reservours
            continue
        if ic.dropna().shape[0] == 0 == 0:
            print(f'{c} has no values in inflow')
            continue
        try:
            # check if country has wind or solar production
            pc = prod_total.sel(country=c).to_dataframe()
            df = pd.concat([dc, ic, pc], axis=1)
            df = df.drop(['country', 'temp'], axis=1)
            ## change demand unit GWh to MWh
            df['demand'] = df['demand']*1e3
            # get demand minus wind and solar production
            # cut of at 0 (no storage)
            df['netto_demand'] = (df['demand']-df['wind_solar_ror'])
            df.loc[df['netto_demand']<0, 'netto_demand'] = 0
        except:
            print(f'{c} has no wind or solar production')
            # if not, just take demand
            df = pd.concat([dc, ic])
            df = df.drop(['country', 'temp'], axis=1)
            ## change demand unit GWh to MWh
            df['demand'] = df['demand']*1e3
            df = df.rename(columns = {'demand': 'netto_demand'}) 
        caph = df_cap.loc[c][0]
        w = best_storage.loc[c, 'week']
        ff = fill_frac.loc[c, fill_frac.columns[0]]
        model, outflow = rm.optimize_reservoirs_usage(df, caph, stor_weeks=w, begin_fraction = ff)
        countries.append(c)
        appended_data.append(outflow)

    res_outflow = pd.concat(appended_data, 
                        keys=countries, 
                        names=['country', 'time'])
    ### save outflow
    res_outflow = res_outflow.to_xarray()
    res_outflow['country'] = [mapping_back[c] for c in res_outflow.country.values]
    ofile = dir_res_out + os.path.basename(demand_file).replace('demand','outflow')
   
    res_outflow = attributes.set_global_attributes(
        res_outflow, 'hydro inflow, demand, PV production, windproduction')
    res_outflow = attributes.set_hydro_attributes(res_outflow)
    res_outflow.time.attrs.update(demand.time.attrs)
    res_outflow = res_outflow.transpose('time', 'country')
    res_outflow.to_netcdf(ofile)
