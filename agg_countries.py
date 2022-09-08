import xarray as xr
import config
import glob
import os
import sys
import numpy as np
# set function path
sys.path.append(config.dir_functions)

import data_processing.masking as masking
import constants.mappings as mappings

# open mappings 
df_countries = mappings.df_countries_select
# demand uses EU_map so drop nan values
df_countries = df_countries.loc[df_countries.EU_map.notna()]
shapefile_countries = config.shapefile_countries


def sum_per_region(region_mask, da):
    """
    sum variables over region

    parameters
    ----------
    region_mask (xarray.Dataset): netcdf region mask
    ds (xarray.Dataarray):  dataarray to sum over regions

    returns
    -------
    dataset (xarray.DataSet) : dataSet with variables summed over region
    """
    dfc = (da * region_mask.mask).sum(dim=["lat", "lon"])
    dfc.attrs = da.attrs
    return dfc

config.print_header("AGGREGATING COUNTRIES")

# get runs
file_extensions = config.runs
# input dirs
dir_pv_wind_in = config.dir_out + "energy_prod_wind_solar/"
dir_res_in = config.dir_out + "energy_inflow_reservoirs/"
dir_ror_in = config.dir_out + "energy_prod_ror/"
# output dirs
dir_agg = config.dir_out + "agg_production/"
dir_pv_wind_agg = dir_agg + "energy_prod_wind_solar/"
dir_res_agg = dir_agg + "energy_inflow_reservoirs/"
dir_ror_agg = dir_agg + "energy_prod_ror/"
# make output folders if they don't exist yet
for d in [dir_agg, dir_pv_wind_agg, dir_res_agg, dir_ror_agg]:
    isExist = os.path.exists(d)
    if not isExist:
        # create a new directory when it does not exist
        os.makedirs(d)

# loop over different production types 
for ofolder, ifolder in zip(
    [dir_pv_wind_agg , dir_res_agg, dir_ror_agg], 
    [dir_pv_wind_in, dir_res_in, dir_ror_in]):
    # open all files in input folder
    files = glob.glob(f'{ifolder}*nc')
    files.sort()
    # use first file to make region mask
    region_mask = masking.make_region_mask(
        df_countries.EU_map,
        files[0],
        shapefile_countries,
        country_indexes=df_countries.index_nr,
    )
    
    # make a list with only the files of the selected runs in it
    files2 = []
    for r in config.runs:
        file = [f for f in files if r in f][0]
        files2.append(file)
    
    # loop over those files
    for f in files2:     
        print(f)
        # open file
        dst = xr.open_dataset(f)
        # loop over years in the file to aggregate
        dst_yearlist = []
        for y in np.unique(dst['time.year'].values):
            print(y)
            ds = xr.Dataset()
            dst2 = dst.sel(time=str(y))
            # loop over variables
            for v in dst2:
            # get mean per region
                dsc = sum_per_region(region_mask, dst2[v])
                print(v, ' is done')
                ds[v] = dsc
                dsc.close()
            dst_yearlist.append(ds)
        # concat years
        ds = xr.concat(dst_yearlist, dim='time')    
        ofile = os.path.basename(f)[:-3] + "_agg.nc"
    
        ds.to_netcdf(
            ofolder+ ofile)
        print(ofolder+ ofile)
        dst.close()
        dsc.close()
    print('grouping is done')
print('run is done')

