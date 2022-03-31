import xarray as xr
import config
import glob
import os
import sys
# set function path
sys.path.append(config.dir_functions)

import data_processing.masking as masking  
import constants.mappings as mappings

def sum_per_region(region_mask, ds):

    dfc = (ds * region_mask.mask ).sum(dim=['lat', 'lon'])
    dfc.attrs = ds.attrs
    return dfc

def Filter(string, substr):
    'mask filters in strings with a list of substrings'
    return [str for str in string if
             any(sub in str for sub in substr)]

config.print_header('AGGREGATING COUNTRIES')

file_extensions = config.runs
# input dirs
dir_pv_wind_in = config.dir_out + 'energy_prod_wind_solar/'
dir_res_in = config.dir_out + 'energy_inflow_reservoirs/'
dir_ror_in = config.dir_out + 'energy_prod_ror/'
# output dirs
dir_agg = config.dir_out + 'agg_production/'
dir_pv_wind_agg = dir_agg + 'energy_prod_wind_solar/'
dir_res_agg = dir_agg + 'energy_inflow_reservoirs/'
dir_ror_agg = dir_agg + 'energy_prod_ror/'

for d in [dir_agg, dir_pv_wind_agg, dir_res_agg, dir_ror_agg]:
    isExist = os.path.exists(d)
    if not isExist:
        # create a new directory when it does not exist
        os.makedirs(d)
res_files = Filter(glob.glob(dir_res_in+"*.nc"), file_extensions)
prod_files = Filter(glob.glob(dir_pv_wind_in+"*.nc"), file_extensions)
ror_files = glob.glob(dir_ror_in+"*.nc")

df_countries = mappings.df_countries_select
# demand uses EU_map so drop nan values 
df_countries = df_countries.loc[df_countries.EU_map.notna()]
shapefile_countries = config.shapefile_countries

region_mask_ror = masking.make_region_mask(
    df_countries.EU_map, 
    ror_files[0], 
    shapefile_countries,
    country_indexes=df_countries.index_nr)    
    
for f in ror_files:
    print('Aggregating ror production:', os.path.basename(f)[:-3])
    ds = xr.open_dataset(f)
    ds = sum_per_region(region_mask_ror, ds)
    ofile = os.path.basename(f)[:-3] + '_agg.nc'
    ds.to_netcdf(dir_ror_agg + ofile)
    ds.close()
    
region_mask_res = masking.make_region_mask(
    df_countries.EU_map, 
    res_files[0], 
    shapefile_countries,
    country_indexes=df_countries.index_nr)

for f in res_files:
    print('Aggregating reservoir inflow:', os.path.basename(f)[:-3])
    ds = xr.open_dataset(f)
    ds = sum_per_region(region_mask_res, ds)
    ofile = os.path.basename(f)[:-3] + '_agg.nc'
    ds.to_netcdf(dir_res_agg + ofile)
    ds.close()
    
region_mask_prod = masking.make_region_mask(
    df_countries.EU_map, 
    prod_files[0], 
    shapefile_countries,
    country_indexes=df_countries.index_nr)    
    
for f in prod_files:
    print('Aggregating PV & wind production:', os.path.basename(f)[:-3])
    ds = xr.open_dataset(f)
    ds = sum_per_region(region_mask_prod, ds)
    ofile = os.path.basename(f)[:-3] + '_agg.nc'
    ds.to_netcdf(dir_pv_wind_agg + ofile)
    ds.close()


