import xarray as xr
import numpy as np
import geopandas as gpd
import regionmask
import json

import os

import data_processing.attributes as attributes
import constants.ldd as c_ldd

import config

area = dict(lon=slice(config.l, config.r), lat=slice(config.t, config.b))


def cut_box(dataset, lon_360=True):
    """
    cuts dataset to european box for gaussion 80 grid
    for lon_360 'True' the longitude is first converted to -180 to 180 degrees

    parameters
    ----------
    dataset (xarrayDataset): dataset to cut
        if 'False' nothing

    returns
    -------
    mask (xarrayDataArray): dataset cut to europe box with -180 to 180 degrees longitude
    """
    if lon_360:
        dataset = dataset.assign_coords(lon=(((dataset.lon + 180) % 360) - 180)).sortby(
            "lon"
        )

    if dataset.lat[0] < dataset.lat[-1]:
        dataset = dataset.loc[
            dict(lon=slice(config.l, config.r), lat=slice(config.b, config.t))
        ]
    else:
        dataset = dataset.loc[
            dict(lon=slice(config.l, config.r), lat=slice(config.t, config.b))
        ]
    return dataset


def _make_region_mask(
    geodf,
    dataset,
    region_selection,
    country_name_column,
    abbreviation_column,
    country_indexes,
    no_numbers=True,
):
    """
    takes geodataframe and country selection to mask dataset from http://www.naturalearthdata.com/downloads/10m-cultural-vectors/
       no_numbers is True: a single mask for the regions
        with value '1' for mask and 'nan' for out of region.
    no_numbers is False: a mask with index per region.
        Can be used to mask different countries.

    parameters
    ----------
    geodf (GeoDataFrame): geodataframe with multipolygons of countries
    dataset (xarrayDataset): dataset to be masked to list of countries
    region_selection (list of strings): list with abbriviations of countries
        to include in masking

    default parameters
    ----------
    country_name_column (string): name of the column that contains
        the country names
    abbreviation_column (string): name of the column that contains
        the abbriviation of the country names
    no_numbers (boolean): 'True' a mask without country indexes,
        'False' a mask with seperate country indexes

    returns
    -------
    mask (xarrayDataArray): mask with 'nan' values out of region and integers for regions
    """

    # select polygons from geodf
    geodf_selec = geodf.copy().loc[geodf[abbreviation_column].isin(region_selection)]

    # if the index_nr are note defined that geodf indec numbers
    if country_indexes is None:
        geodf_selec["index_nr"] = list(geodf_selec.index)
    else:
        geodf_selec["index_nr"] = list(
            map(
                dict(zip(region_selection.values, country_indexes.values)).get,
                geodf_selec[abbreviation_column],
            )
        )

    # make polygon masks
    countries_mask_poly = regionmask.Regions(
        name="countries",
        numbers=geodf_selec.index_nr,  # indexes,
        names=geodf_selec[country_name_column],  # geodf[country_name_column][indexes],
        abbrevs=geodf_selec[abbreviation_column],
        outlines=list(
            x for x in geodf_selec.geometry.values
        ),  # list(geodf.geometry.values[i] for i in indexes)
    )

    if "time" in dataset.dims:
        layer = dataset.isel(time=0)
    else:
        layer = dataset

    # make dataset mask based on polys
    mask = countries_mask_poly.mask(layer, lat_name="lat", lon_name="lon")

    if no_numbers:
        mask = xr.where(mask >= 0, 1, mask)

    # update variable attributes
    mask.attrs.update(
        standard_name="region mask",
        long_name="Boolean mask for selected region(s)",
        units="[-]",
    )

    # update dimensions attributes
    #     try:
    #         mask = attributes.set_lat_lon_attributes(mask, list(layer.dims.keys()))
    #     except:

    #         for x in layer.dims:
    #             mask[x].attrs.update(layer[x].attrs)

    return xr.Dataset({"mask": mask})


def _split_mask_into_regions(mask):
    """
    splits a mask with dimensions (lat,lon) and different integers per region in a mask with 'nan'
        values out of region and int(1) in region, and the dimenions (lat,lon,country)

    parameters
    ----------
    mask (xarrayDataArray): mask with 'nan' values out of region and integers for regions.
        With dimensions (lat,lon)

    returns
    -------
      splitted_mask (xarrayDataArray): mask with 'nan' values out of region and int(1) for regions.
        With dimensions (lat,lon, region)
    """
    xr.set_options(keep_attrs=True)
    mask = mask.mask
    splitted_mask = xr.concat(
        [mask.where(mask == c) for c in np.unique(mask.data[~np.isnan(mask.data)])],
        "country",
    )
    splitted_mask["country"] = np.unique(mask.data[~np.isnan(mask.data)])
    # remove masking numbers
    splitted_mask = xr.where(splitted_mask.notnull(), 1, np.nan)
    splitted_mask = xr.Dataset({"mask": splitted_mask})
    return splitted_mask


def _make_inflows_mask_from_mask(region_mask, ldd, change=c_ldd.pcr_ldd_directions):
    """
    Based on a region_mask with with 'nan' values out of region and int(1) for regions and
        dimensions (lat,lon,region) make a mask of a 1-pixel wide outline of the region mask and
        select the cells with a flow direction into the country

    parameters
    ----------
    region_mask (xarrayDataArray): mask with 'nan' values out of region and integers for regions.
        with dimensions (lat,lon, region)
    ldd (xarrayDataSet): local drain direction map directions

    default parameters
    ----------
    change (dict): dictionary of flow direction in steps [x,y] per int in ldd. With str(int) as keys.

    returns
    -------
      (xarrayDataSet): mask with 'nan' values out of inflow and int(1) for cells where drain
      direction enters the region
    """
    # select dataArray
    region_mask = region_mask.mask
    # make inverse mask and use to mask ldd
    mask_inverse = xr.where(region_mask == 1, 0, 1)
    border_flow = (mask_inverse * ldd).where((mask_inverse * ldd) > 0)

    # set empty dataArray
    ds_empty = xr.where(border_flow == 100, 1, np.nan)
    lons = ds_empty.lon
    lats = ds_empty.lat

    dst2 = ds_empty

    # loop over different drain directions, and move from according to drain direction
    for i in [int(x) for x in change.keys()]:
        # move flow at border in 1 direction
        dst1 = xr.where(border_flow == i, i, np.nan)
        dst1["lon"] = lons + (0.5 * change[str(i)][0])
        dst1["lat"] = lats + (0.5 * change[str(i)][1])
        #     Only keep the ldd if it enters the region
        dst1 = dst1.where(region_mask > 0)
        #     move the ldd-cell back to the border
        dst1["lon"] = dst1.lon.data - (0.5 * change[str(i)][0])
        dst1["lat"] = dst1.lat.data - (0.5 * change[str(i)][1])
        #   concat all drain directions
        dst2 = xr.concat([dst2, dst1], "dir")
    # sum all drain directions
    inflows = dst2.sum(dim="dir")
    inflows = inflows.rename({"ldd": "mask"})

    # update dimensions attributes
    try:
        inflows = attributes.set_lat_lon_attributes(inflows, list(ldd.dims.keys()))
    except:
        for x in ldd.dims:
            inflows[x].attrs.update(ldd[x].attrs)

    return xr.where(inflows, 1, np.nan)


def make_region_mask(
    region_selection,
    gridfile,
    shapefile_countries,
    country_name_column=config.country_name_column,  # config'ADMIN',
    abbreviation_column=config.abbreviation_column,  #'ADM0_A3',
    # europe=True,
    country_indexes=None,
):

    """
    masks dataset based on a list of regions and a region shapefile
    geodataframe and country selection to mask dataset from:
    http://www.naturalearthdata.com/downloads/10m-cultural-vectors/

    parameters
    ----------
    region_selection (list of strings): list with abbriviations of countries
        to include in masking
    gridfile (xarrayDataset): dataset with grid to be masked to list of countries
    ofile (xarrayDataArray): dataset with new dimension of regions that have been masked
    shapefile_countries (shp): shapefile that contains the regions that need to be masked

    default parameters
    ----------
    country_name_column (string): name of the column that contains
        the country names
    abbreviation_column (string): name of the column that contains
        the abbriviation of the country names
    europe (boolean): 'True' cut europe to reduce processing time,
        'False' full dataset
    country_indexes(list/PandaSeries of int): list to number countries

    returns
    -------
    dsmasked (xarrayDataArray): mask
    """
    xr.set_options(keep_attrs=True)
    # open datasets
    regions = gpd.read_file(shapefile_countries)

    # if marine eex areas only select eez and landlocked countries, ignore double regions
    try:
        regions = regions.loc[
            regions.POL_TYPE.isin(["Union EEZ and country", "Landlocked country"])
        ]
    #         regions = regions.loc[regions['POL_TYPE'].isin(['Union EEZ and country','Landlocked country'])]
    except:
        pass

    ds = xr.open_dataset(gridfile)
    # rename lon and latitude and cut to box to increase processing time
    try:
        ds = cut_box(ds)
    except:
        ds = ds.rename({"longitude": "lon", "latitude": "lat"})
        ds = cut_box(ds)

    # mask to dataset per countries
    mask = _make_region_mask(
        regions,
        ds,
        region_selection,
        country_name_column,
        abbreviation_column,
        country_indexes=country_indexes,
        no_numbers=False,
    )

    mask = _split_mask_into_regions(mask)

    return mask


def cut_netcdf_into_regions(
    region_selection,
    ifile,
    ofile,
    shapefile_countries,
    country_name_column=config.country_name_column,  # config'ADMIN',
    abbreviation_column=config.abbreviation_column,  #'ADM0_A3',
    # europe=True,
    country_indexes=None,
):

    """
    masks dataset based on a list of regions and a region shapefile
    geodataframe and country selection to mask dataset from:
    http://www.naturalearthdata.com/downloads/10m-cultural-vectors/

    parameters
    ----------
    region_selection (list of strings): list with abbriviations of countries
        to include in masking
    ifile (xarrayDataset): dataset to be masked to list of countries
    ofile (xarrayDataArray): dataset with new dimension of regions that have been masked
    shapefile_countries (shp): shapefile that contains the regions that need to be masked

    default parameters
    ----------
    country_name_column (string): name of the column that contains
        the country names
    abbreviation_column (string): name of the column that contains
        the abbriviation of the country names
    europe (boolean): 'True' cut europe to reduce processing time,
        'False' full dataset
    country_indexes(list/PandaSeries of int): list to number countries

    returns
    -------
    dsmasked (xarrayDataArray): masked dataset
    """

    mask = make_region_mask(
        region_selection,
        ifile,
        shapefile_countries,
        country_name_column=config.country_name_column,  # config'ADMIN',
        abbreviation_column=config.abbreviation_column,  #'ADM0_A3',
        # europe=True,
        country_indexes=country_indexes,
    )

    ds = xr.open_dataset(ifile)
    try:
        ds = cut_box(ds)
    except:
        ds = ds.rename({"longitude": "lon", "latitude": "lat"})
        ds = cut_box(ds)
    ds_countries = mask.mask * ds

    ds_countries.to_netcdf(ofile, "w")

    if country_indexes is None:
        # also add country names to integers
        regions = gpd.read_file(shapefile_countries)
        country_list = list(regions[abbreviation_column])
        indexes = [country_list.index(x) for x in region_selection]
        indexes.sort()
        ordered_countries = [regions[abbreviation_column].iloc[x] for x in indexes]
        # create json object from remap dictionary
        map_of_names = dict(zip(indexes, ordered_countries))
        fname = (
            os.path.dirname(ofile) + "/dict_" + os.path.basename(ofile)[:-3] + ".json"
        )

        with open(fname, "w") as fp:
            json.dump(map_of_names, fp)

    else:
        map_of_names = dict(zip(country_indexes, region_selection))
        fname = (
            os.path.dirname(ofile) + "/dict_" + os.path.basename(ofile)[:-3] + ".json"
        )

        with open(fname, "w") as fp:
            json.dump(map_of_names, fp)


# def make_inflow_mask(region_selection,
#                     lddfile,
#                     shapefile_countries,
#                     country_name_column = config.country_name_column, #config'ADMIN',
#                     abbreviation_column = config.abbreviation_column, #'ADM0_A3',
#                     europe=True,
#                     country_indexes = None):
#     """
#     Based on a region_mask with with 'nan' values out of region and int(1) for regions and
#         dimensions (lat,lon,region) make a mask of a 1-pixel wide outline of the region mask and
#         select the cells with a flow direction into the country

#     parameters
#     ----------
#     region_mask (xarrayDataArray): mask with 'nan' values out of region and integers for regions.
#         with dimensions (lat,lon, region)
#     ldd_file (str): directory of netcdf file with local drain directions

#     default parameters
#     ----------
#     change (dict): dictionary of flow direction in steps [x,y] per int in ldd. With str(int) as keys.

#     returns
#     -------
#       (xarrayDataArray): mask with 'nan' values out of inflow and int(1) for cells where drain
#       direction enters the region
#     """
#     xr.set_options(keep_attrs=True)
#     # open datasets
#     regions = gpd.read_file(shapefile_countries)
#     ds = xr.open_dataset(lddfile)

#     # rename lon and latitude and cut to box to increase processing time
#     try:
#         ds = cut_box(ds)
#     except:
#         ds = ds.rename({'longitude':'lon','latitude':'lat' })
#         ds = cut_box(ds)

#     mask = _make_region_mask(regions,
#                              ds,
#                              region_selection,
#                              country_name_column,
#                              abbreviation_column,
#                              country_indexes=country_indexes,
#                              no_numbers = False,
#                             )

#     region_mask = _split_mask_into_regions(mask)
#     inflow_mask = _make_inflows_mask_from_mask(region_mask,ds)
#     return inflow_mask


# def mask_countries(dataset, mask, indexes):
#     """
#     masks dataset based on list of integers when multiple countries in mask

#     parameters
#     ----------

#     dataset (xarrayDataset): dataset to be masked to list of countries
#     mask (xarrayDataArray): mask with 'nan' values out of region and integers for regions
#     indexes (list of int): list of integers to select regions to mask

#     returns
#     -------
#     dsmasked (xarrayDataArray): masked dataset
#     """
#     for country in indexes:
#         temp = dataset.where(mask==country)
#         if country == indexes[0]:
#             dsmasked = temp
#         else:
#             dsmasked = xr.merge([dsmasked, temp])

#     return dsmasked
