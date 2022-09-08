import os
import glob
import config
import sys
import time 

import numpy as np
import pandas as pd
import xarray as xr
import pyomo.environ as pyo

# set function path and import functions
sys.path.append(config.dir_functions)
import constants.mappings as mappings
import optimize_reservoirs.reservoir_model_rolling as rm
import data_processing.attributes as attributes


def pad_list(some_list, target_len):
    """
    make a list the desired length

    Parameters
    ----------
    some_list (list) : 
    target_len(int): 

    Returns
    -------
    list padded to be target length 
    """
    return some_list[:target_len] + [0]*(target_len - len(some_list))

## country mapping
df_countries = mappings.df_countries_select

### imort variables
reservoir_size = pd.read_csv(config.storage_capacities, index_col=0)
reservoir_size = reservoir_size.replace(np.nan, 1)
fill_frac = pd.read_csv(config.filling_level_at_start, index_col=0)

## Run optimization settings
# set overlap between windows and size of windows
overlap =14
window_size = 28
stepsize=window_size+overlap

# Create a solver
opt = pyo.SolverFactory("glpk") #('ipopt')

### open files
## mappings
df_countries = mappings.df_countries_select
# demand uses EU_map so drop nan values
df_countries = df_countries.loc[df_countries.EU_map.notna()]
# make mappings
mapping = dict(zip(df_countries.index_nr, df_countries.EU_map))
mapping_cap = dict(zip(df_countries.nuts_id, df_countries.index))
mapping_back = dict(zip(df_countries.EU_map, df_countries.index_nr))
## input directories
dir_bounds = config.dir_out + "proccessed_bounds/"
dir_agg = config.dir_out + "agg_production/"
dir_demand = config.dir_out + "energy_demand/"
dir_pv_wind = dir_agg + "energy_prod_wind_solar/"
dir_res_in = dir_agg + "energy_inflow_reservoirs/"
dir_ror = dir_agg + "energy_prod_ror/"

### set output directories
dir_res_out = dir_agg + "energy_outflow_reservoirs_rolling/"
# make output directory
for d in [dir_res_out]:
    isExist = os.path.exists(d)
if not isExist:
    # create a new directory when it does not exist
    os.makedirs(d)

# open country capacities based on gridded input
reservoir_cap = config.hydro_cap
df_cap = xr.open_dataset(reservoir_cap).to_dataframe()
df_cap = df_cap.groupby("country_code").sum()
df_cap = df_cap[["capacity"]]
# clean data
df_cap.replace("", np.nan, inplace=True)
df_cap.replace(0, np.nan, inplace=True)
df_cap.dropna(inplace=True)
# map country_code
df_cap.index = [mapping_cap[k] for k in df_cap.index]

### open files
demand_files = glob.glob(dir_demand + "*.nc")
inflow_files = glob.glob(dir_res_in + "*.nc")
production_files = glob.glob(dir_pv_wind + "*.nc")
ror_prod_files = glob.glob(dir_ror + "*.nc")

### start optimmization
config.print_header(f"COMPUTING HYDRO OUTFLOW")
# loop over runs in ensemble
for r in config.runs:
    print("Computing dispatch for", r)
    # set infiles
    demand_file = [f for f in demand_files if r in f][0]
    inflow_file = [f for f in inflow_files if r in f][0]
    production_file = [f for f in production_files if r in f][0]
    ror_file = [f for f in ror_prod_files if r in f][0]

    # demand
    demand = xr.open_dataset(demand_file)
    demand["country"] = list(map(mapping.get, demand["country"].values))
    demand = demand.dropna(dim="time")
    # inflow
    inflow = xr.open_dataset(inflow_file)
    inflow["country"] = list(map(mapping.get, inflow["country"].values))
    # production
    prod_m = xr.open_dataset(production_file)
    prod_m["country"] = list(map(mapping.get, prod_m["country"].values))
    # run of river production
    ror_prod = xr.open_dataset(ror_file)
    ror_prod["country"] = list(map(mapping.get, ror_prod["country"].values))
    # sum all vars in production (wind and solar)
    ror_prod = ror_prod.loc[dict(time=ror_prod.time.where(demand.time))]
    prod_m["ror"] = ror_prod.ror

    # sum all vars in production (wind, solar and ror)
    vars_to_sum = list(prod_m.keys())
    prod_total = (
        prod_m[vars_to_sum].to_array().sum("variable").to_dataset(name="wind_solar_ror")
    )

    # initiate lists for loop
    appended_data = []
    countries = []
    # loop over all countries for which demand was computes
    for c in demand.country.values:
        print('run:', r, 'running for ', c)
        # open country demand
        dc = demand.sel(country=c).sortby("time").to_dataframe()
        try:
            # open reservoir inflow and drop 0 values
            ic = inflow.sel(country=c).to_dataframe()
            ic = ic.replace(0, np.nan)
        except:
            # no need to optimize if country has no hydro reservours
            print(f"{c} not in inflow file")
            continue
        if ic.dropna().shape[0] == 0 == 0:
            # empy inflow array so no optimization
            print(f"{c} has no values in inflow")
            continue
        try:
            # check if country has wind or solar production
            pc = prod_total.sel(country=c).to_dataframe()
            df = pd.concat([dc, ic, pc], axis=1)
            df = df.drop(["country", "temp"], axis=1)
            ## change demand unit GWh to MWh
            df["demand"] = df["demand"] * 1e3
            # get demand minus wind and solar production
            # cut-off at 0 (no storage)
            df["netto_demand"] = df["demand"] - df["wind_solar_ror"]
            df.loc[df["netto_demand"] < 0, "netto_demand"] = 0
        except:
            print(f"{c} has no wind or solar production")
            # if not, just take demand
            df = pd.concat([dc, ic])
            df = df.drop(["country", "temp"], axis=1)
            ## change demand unit GWh to MWh
            df["demand"] = df["demand"] * 1e3
            df = df.rename(columns={"demand": "netto_demand"})
        # generate mean yearly cycle of input values
        dfmean = df.groupby(df.index.dayofyear).mean()

        ### get country specific boundary conditions
        # get hydropower capacity to bound outflow
        caph = df_cap.loc[c][0]
        # get reservoir storage size for country
        storage_size_in_weeks = reservoir_size.loc[c, "week"]
        stor_size = df.inflow.resample('1Y').sum()[1:].mean()/(52/storage_size_in_weeks) # convert to MWh
        # get inital fill fraction of country reservoir 
        ff = fill_frac.loc[c, fill_frac.columns[0]]
        if ff>1:
            ff=1

        # initiate list to append rolling windows 
        windows = []


        ### get estimated yearly cycle for output variables
        # run optimization once with perfec foresight
        model, outflow0, results = rm.optimize_reservoirs_usage(
                    df, caph, stor_size, begin_fraction=ff
                )
        if (results.solver.termination_condition != 'optimal'):
                print(i,'!!!!--NO FEASIBLE SOLUTION--!')
        # get mean yearly reservoir cycle from result   
        rescycle = outflow0.reservoir.groupby(outflow0.index.dayofyear).mean()

        # get intital timestamp and fillfactor to start first window
        ts = outflow0.iloc[0].name
        ff = outflow0.iloc[0].reservoir/stor_size

        # make a 1 year prediction of mean input variables 
        long_term_timerange = pd.date_range(start=ts, end=ts+pd.offsets.DateOffset(years=1), freq='1D')
        dflong_term = dfmean.loc[long_term_timerange.dayofyear,:]
        dflong_term.index = long_term_timerange
        # optimize for low inflow, so storage will be filled.
        dflong_term['inflow'] = dflong_term.inflow*0.5
        # replace the first window with short term predictions
        ff2 = rescycle.loc[dflong_term.index[-1].dayofyear]/stor_size
        dflong_term.iloc[0:stepsize, :] = df.loc[dflong_term.index[0:stepsize]]

        # run the first rolling optimzation
        model, outflow1, results = rm.optimize_reservoirs_usage(
                    dflong_term, caph,stor_size, begin_fraction=ff, end_fraction=ff2
                )

        # append the outflow to the windows results
        windows.append(outflow1)

        # loop optimization over rolling windows for as long as the run is
        for i in range(1,(len(df)//(stepsize-overlap)-1)):
            ### get initial conditions
            # timestep of outflow where next window begins
            ts = outflow1.iloc[stepsize-overlap].name
            # fill fraction of outflow where next window begins
            ff = outflow1.iloc[stepsize-overlap].reservoir/stor_size
            # fill fraction at the end of previous window (overlap)
            ffbound = outflow1.iloc[stepsize].reservoir/stor_size
            startflow = float(outflow1.iloc[stepsize-overlap].Eout)
 
            # make a 1 year prediction of mean input variables starting at ts
            long_term_timerange = pd.date_range(start=ts, end=ts+pd.offsets.DateOffset(years=1), freq='1D')
            dflong_term = dfmean.loc[long_term_timerange.dayofyear,:]
            dflong_term.index = long_term_timerange
            # replace the first window with short term predictions
            dflong_term.iloc[0:stepsize, :] = df.loc[dflong_term.index[0:stepsize]]
            # set mean fill fraction at end of run as boundary condition
            ff2 = rescycle.loc[dflong_term.index[-1].dayofyear]/stor_size

            # run the rolling optimization with boundary conditions
            model, outflow1, results = rm.optimize_reservoirs_usage(
                        dflong_term, caph, stor_size, begin_fraction=ff, end_fraction = ff2, 
                init=pad_list(list(outflow1.Eout.values), len(dflong_term))[0:len(dflong_term)], 
                startflow=startflow, overlap=overlap, bound_fraction = ffbound
                    )
            # check if solution is feasible
            if (results.solver.termination_condition != 'optimal'):
                print(i,'!!!!--NO FEASIBLE SOLUTION--!')
            # append the outflow to the windows results
            windows.append(outflow1)
        # concat the windows of all rolling optimizations in a dataframe
        country_output = pd.concat([w[overlap:stepsize] for w in windows])
        # append the output of each country to a list
        countries.append(c)
        appended_data.append(country_output)
    # concat the country results to a pandas dataframe
    res_outflow = pd.concat(appended_data, keys=countries, names=["country", "time"])
    # export to xarray
    res_outflow = res_outflow.to_xarray()
    res_outflow["country"] = [mapping_back[c] for c in res_outflow.country.values]
    # set attributes
    res_outflow = attributes.set_global_attributes(
        res_outflow, "hydro inflow, demand, PV production, windproduction"
    )
    res_outflow = attributes.set_hydro_attributes(res_outflow)
    res_outflow.time.attrs.update(demand.time.attrs)
    res_outflow = res_outflow.transpose("time", "country")
    # save output file
    ofile = dir_res_out + os.path.basename(demand_file).replace('demand','outflow')[:-3] + '_rolling.nc'
    res_outflow.to_netcdf(ofile)
    print(f'run {r} is done')
print('hydropower dispatch is done')