import xarray as xr
import pcraster as pcr
import numpy as np
from cdo import Cdo

cdo = Cdo()  # python
import data_processing.masking as masking
import constants.ldd as c_ldd
import data_processing.attributes as at

steps = c_ldd.pcr_ldd_directions


def define_outflows_of_cutout(localDrainDirectionMap):
    """
    replace all cells that flow out of the map area with flow direction '5'.
        This function only works with pcr notation of ldd directions.

    parameters
    ----------
    localDrainDirectionMap (xarray.Dataset): data_set with local drain directions

    returns
    -------
    (xarray.Dataset): data array for the land area of europe with outlfow cells with
        flow directionn '5'
    """

    ds = localDrainDirectionMap
    ## replace all outflows at end of map with '5' (no flow direction)
    t1 = xr.where(
        (ds.lat == ds.lat.min().data) & ((ds.ldd == 3) | (ds.ldd == 2) | (ds.ldd == 1)),
        5,
        ds,
    )

    t1 = xr.where(
        (ds.lat == ds.lat.max().data) & ((ds.ldd == 7) | (ds.ldd == 8) | (ds.ldd == 9)),
        5,
        t1,
    )

    t1 = xr.where(
        (ds.lon == ds.lon.min().data) & ((ds.ldd == 1) | (ds.ldd == 4) | (ds.ldd == 7)),
        5,
        t1,
    )

    t1 = xr.where(
        (ds.lon == ds.lon.max().data) & ((ds.ldd == 9) | (ds.ldd == 6) | (ds.ldd == 3)),
        5,
        t1,
    )

    return t1.transpose("lat", "lon")


def clone_attributes():
    """
    Get the map attributes from the PCRaster clone map.

    Returns a list with: xmin, xmax, ymin, ymax, nr_rows, nr_cols, cell_size.
    """
    nr_rows = pcr.clone().nrRows()
    nr_cols = pcr.clone().nrCols()
    cell_size = pcr.clone().cellSize()
    ymax = pcr.clone().north()
    xmin = pcr.clone().west()

    return xmin, ymax, nr_rows, nr_cols, cell_size


def pcr_2D_to_netcdf(
    map, mappcr, varname, attributes=None, no_negative=True, unit="[-]"
):
    """
    read pcr.map file in xarray dataset

    parameters
    ----------
    map (str): filename of 2D raster in pcraster format to set as clone (xxxx.map)
    mappcr (ppcraster._pcraster.Field): rastermap with data to save in netcdf
    varname (str): name of data for var in xarray.Dataset
    attributes (dict): global attributes to add to dataset


    returns
    -------
    (xarray.Dataset): dataset with raster data as one of the variables
    """

    pcr.setclone(map)
    # read coordiantes of pcr and make lat and lon arrays
    xmin, ymax, nr_rows, nr_cols, cell_size = clone_attributes()
    lons = np.arange(xmin + 0.5 * cell_size, xmin + (cell_size * nr_cols), cell_size)
    lats = np.arange(ymax - (0.5 * cell_size), ymax - (cell_size * nr_rows), -cell_size)

    # conver pcr to numpy and replace 'misisng values' with nan
    if no_negative:
        MV = 1e21
        data = pcr.pcr2numpy(mappcr, MV)
        data = np.where(data > 0, data, np.nan)
        data = np.where(data != MV, data, np.nan)
    else:
        MV = -9999
        data = pcr.pcr2numpy(mappcr, MV)
        data = np.where(data != MV, data, np.nan)

    # make dataset
    ds = xr.Dataset(
        data_vars={varname: (["lat", "lon"], data)},
        coords=dict(
            lon=(["lon"], lons),
            lat=(["lat"], lats),
        ),
        attrs=attributes,
    )

    try:
        ds.attrs.update(attributes)
    except:
        ds = at.set_global_attributes(ds, "pcraster ldd", grid="", area="")

    # update variables and coordinate attributes
    ds = at.set_lat_lon_attributes(ds)

    ds[varname] = ds[varname].assign_attrs(
        dict(
            standard_name=varname,
            units=unit,
        )
    )
    return ds


def make_routing_layer_dataset_from_ldd(ds_ldd):
    """
    convert the local drain direction map in a dataset with layers for each routing step.

    parameters
    ----------
    ds_ldd (xarray.DataSet): 2D dataset (lat,lon) with local drain direction

    returns
    -------
    (xarray.Dataset): a dataset with the dimensions lat, lon and layer with layer=0 containing the
        cell furtherst from the sea and layer=-1 containing all the outflows.
    """

    # make an empty dataset
    ds_empty = xr.where(ds_ldd == 100, np.nan, np.nan)
    lons = ds_empty.lon
    lats = ds_empty.lat

    start = ds_ldd.where(ds_ldd == 5)
    ds_layer = start

    n = 150
    for j in range(0, n):
        dst2 = ds_empty
        for i in range(1, 10):
            lats = start.lat.data
            lons = start.lon.data
            dst1 = xr.where(start > 0, i, 0)
            dst1["lon"] = lons + (-0.5 * steps[str(i)][0])
            dst1["lat"] = lats + (-0.5 * steps[str(i)][1])

            dst1 = dst1.where(dst1 == ds_ldd)
            dst1 = dst1.where(dst1 != 5)
            dst2 = xr.concat([dst2, dst1], "dir")

        ds2 = dst2.sum(dim="dir")
        start = ds2
        layer = ds2.where((ds2.ldd > 0))
        ds_layer = xr.concat([ds_layer, layer], "layer")
        if layer.ldd.sum() == 0:
            break
    ds_layer["layer"] = list(reversed(range(0, j + 2)))
    return ds_layer


def _shift_dataset(summed_runoff, da_ldd, direction):
    """
    When routing runoff shift the runoff layer depending on the direction of the ldd

    """
    dst = summed_runoff.where(da_ldd == direction, drop=True)
    dst["lon"] = dst["lon"] + (0.5 * steps[str(direction)][0])
    dst["lat"] = dst["lat"] + (0.5 * steps[str(direction)][1])

    return dst


def _runoff_in_m3_from_file(inputFile, varname):
    """
    Multiply runoff and area of gridcell to m3 for routing
    """
    ds = xr.open_dataset(inputFile)
    # cut to area
    ds = masking.cut_box(ds)
    area = xr.open_dataset(cdo.gridarea(input=inputFile))
    ds["cell_area"] = area.cell_area
    ds["m3"] = ds[varname] * ds.cell_area

    ds["m3"] = ds["m3"].assign_attrs(
        dict(
            standard_name="runoff",
            long_name="total runoff in m3 per gridcell",
            units="m3",
        )
    )

    # drop all other variables
    ds = ds.drop([x for x in ds if x not in ["m3", varname]])

    return ds


def accumulate_runoff(inputFile, ds_ldd, ds_routing, outputFile, varname="mrros"):
    runoff = _runoff_in_m3_from_file(inputFile, varname)
    ds = runoff.where(ds_ldd.ldd > 0).sortby("lat")

    for j in range(0, 70):
        layer = ds_routing.sel(layer=j)
        layer = layer.where(layer > 0, drop=True)
        ds0 = xr.concat(
            [_shift_dataset(ds, layer.ldd, i) for i in range(1, 10)], "dir"
        ).sum(dim="dir")
        ds = xr.concat([ds, ds0], "step").sum(dim="step")
    ds = ds.drop("layer")

    ds = at.set_global_attributes(ds, inputFile, grid="", area="")

    ds = ds.rename({"m3": "discharge"})
    ds["discharge"] = ds["discharge"].assign_attrs(
        dict(
            standard_name="discharge",
            long_name="discharge",
            units="m3/timestep",
        )
    )

    ds = ds.where(ds > 0)

    ds.to_netcdf(outputFile)