# import packages
import xarray as xr
import os
import glob
import time
import sys

## import own modules
import config

# set function path
sys.path.append(config.dir_functions)
import data_processing.masking as masking
import routing.routing as routing

config.print_header("ROUTING RUNOFF")

# ----------------------------------------------------------
### Set files and directories
# ----------------------------------------------------------
dir_out = config.dir_out
dir_bounds = dir_out + "proccessed_bounds/"
dir_disout = dir_out + "ro_accu/"
dir_routout = config.dis_dir
# input files
ldd_file = config.ldd_file
files_ro = config.ro_files
# make output directories
for d in [dir_bounds, dir_disout, dir_routout]:
    isExist = os.path.exists(d)
    if not isExist:
        # create a new directory when it does not exist
        os.makedirs(d)

actions = config.route_runoff_submodules
area = masking.area

var1 = config.ro_varname
# ----------------------------------------------------------
### Make drain direction files
# ----------------------------------------------------------
try:
    # open
    ds_routing = xr.open_dataset(dir_bounds + "lddRoutingLayers.nc")
    ds_ldd = xr.open_dataset(dir_bounds + "lddRegionCutout.nc")
except:
    if actions["make_lddRoutingLayers"]:
        ds_ldd_glob = xr.open_dataset(ldd_file)
        # make cutout from global
        ds_ldd = ds_ldd_glob.loc[area]
        ds_ldd = ds_ldd.where(ds_ldd > 0)
        # define outflows of area
        ds_ldd = routing.define_outflows_of_cutout(ds_ldd)
        # split in routing steps
        ds_routing = routing.make_routing_layer_dataset_from_ldd(ds_ldd)
        # save files
        ds_ldd.close()
        ds_ldd.to_netcdf(dir_bounds + "lddRegionCutout.nc")
        ds_routing.close()
        ds_routing.to_netcdf(dir_bounds + "lddRoutingLayers.nc")

    # open
    ds_routing = xr.open_dataset(dir_bounds + "lddRoutingLayers.nc")
    ds_ldd = xr.open_dataset(dir_bounds + "lddRegionCutout.nc")

# ----------------------------------------------------------
### run for runoff files
# ----------------------------------------------------------

start0 = time.time()
### open runoff files:
for r in config.runs:
    print("Routing runoff for", r)
    # get file names of run
    f = [fro for fro in config.ro_files if r in fro][0]
    start = time.time()
    ofile = dir_disout + os.path.basename(f).replace(var1, "ro_accu")
    routing.accumulate_runoff(f, ds_ldd, ds_routing, ofile, varname=var1)
    end = time.time()

    # write status to txtfile
    update = [
        ofile,
        "time to process:" + str(round(end - start, 3)),
        "total processing time:" + str(round(end - start0, 3)) + "\n",
    ]
    with open("update.txt", "a") as up:
        up.writelines("\n".join(update))


config.print_header("RECESSION COEFFICIENT")

kx = config.kx

for f in glob.glob(dir_disout + "*.nc"):
    dis_accu = xr.open_mfdataset(f)
    da = dis_accu.drop(var1)

    for i, t in enumerate(da.time):
        if i == 0:
            rout = da.isel(time=i)
            new = rout
            new = new.expand_dims("time")
        else:
            rout = (1 - kx) * da.isel(time=i) + kx * rout
            #         rout = rout.expand_dims('time')
            rout["time"] = [t.values]
            new = xr.concat([new, rout], "time")
    fname = os.path.basename(f).replace("ro_accu", "dis")
    new.to_netcdf(dir_routout + fname)
    print(fname)
