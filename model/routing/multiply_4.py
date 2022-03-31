import glob 
import xarray as xr
import os

dir_in =  '/nobackup/users/most/HiWAVES_eu/present/mrros_d_30min/'
files = glob.glob(dir_in+'*.nc')
dir_out =  '/nobackup/users/most/HiWAVES_eu/present/mrros_d_30min4/'

for f in files:
	print(f)
	ds = xr.open_dataset(f)
	ds = 4 * ds

	fout= dir_out+os.path.basename(f)
	ds.to_netcdf(fout)
	ds.close()