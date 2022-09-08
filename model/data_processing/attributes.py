from datetime import datetime
import config


def set_hydro_attributes(ds):
    a = {
        "netto_demand": [
            "netto demand",
            "demand with PV solar and wind production subtracted",
        ],
        "inflow": ["grid cell inflow", "inflow into gridcell"],
        "reservoir": ["reservoir storage", "energy stored in reservoirs at timestep"],
        "Ein": [
            "inflow into reservoir",
            "inflow into reservoir, considering maximum storage (so overflow)",
        ],
        "Eout": ["hydro energy production", "production of reservoir hydro energy"],
        "conventional": [
            "conventional energy production",
            "production not by renewable sources",
        ],
        "diff_dev": [
            "deviation difference",
            "difference of relative mean deviation between demand and hydro production",
        ],
        
        "event": [
            "no production available",
            "energy that can't be produced by renewable or installed conventional sources",
        ],
    }
    for v in ds:
        ds[v].attrs.update(standard_name=a[v][0], long_name=a[v][1], units="MWh")
    return ds


# routing_attributes = dict(
#         description="description",
#         author= config.author_name,
#         data_source = "global local drain direction map for 0.5x0.5 grid",
#         map_area = "global",
#         computation = "computed with pcraster"
#     )


def set_global_attributes(dataset, source, grid="gaussian n80", area="EU13+2"):
    """
    set global attributes to xarray.dataset

    parameters
    ----------
    dataset (xarray.Dataset): dataset for which to update attributes
    parameter (string): option to set for different variables 'demand', 'etc..'

    returns
    -------
    dataset (xarray.DataSet) : dataSet with update attributes
    """

    dataset.attrs.update(
        author=config.author_name,
        project=config.project_name,
        source=source,
        history=f'Computed {datetime.now().strftime("%d-%b-%Y (%H:%M)")}',
        area=area,
        grid=grid,
    )

    return dataset


def set_lat_lon_attributes(dataset, names=["lat", "lon"]):
    """
    set attributes to lat and lon dimensions of xarray.dataset

    parameters
    ----------
    dataset (xarray.Dataset): dataset for which to update attributes
    parameter (string): option to set for different variables 'demand', 'etc..'

    returns
    -------
    dataset (xarray.DataSet) : dataSet with update attributes
    """

    dataset = dataset.rename({names[0]: "lat", names[1]: "lon"})
    dataset.lat.attrs.update(
        standard_name="latitude", long_name="latitude", units="degrees_north", axis="Y"
    )
    
    dataset.lon.attrs.update(
        standard_name="longitude", long_name="longitude", units="degrees_east", axis="X"
    )
    
    return dataset



def set_demand_attributes(dataset):
    dataset.demand.attrs.update(
        standard_name="demand",
        long_name="demand computed from weighted T",
        units="GWh",
    )
    return dataset

# def set_energy_attributes(dataset, Etype):
#     dataset.demand.attrs.update(
#         standard_name = Etype + ' production',
#         long_name = Etype + ' production',
#         units = 'MWh',
#     )
#     return dataset

# def set_hydro_attributes(ds):
#     attrs = {
#         'netto_demand': ['netto demand', 'demand with PV solar and wind production subtracted'],
#         'reservoir': ['reservoir storage', 'energy stored in reservoirs at timestep'],
#         'Ein': ['inflow into reservoir', 'inflow into reservoir, considering maximum storage (so overflow)'],
#         'Eout': ['hydro energy production', 'production of reservoir hydro energy']
#     }
#     for v in ds:
#         ds[v].attrs.update(
#             standard_name = attrs[v][0],
#             long_name =  attrs[v][1],
#             units = 'MWh'
#         )
#     return ds


def update_energy_attributes(
    energy_production, energy_type, unit, eur=False, potential=False
):
    """
    If `eur` is False (default), generates attribute description for location specific , and returns it.
    If `eur` is True, then instead generates attributes for values summed over Europe, and returns that.

    parameters
    ----------
    energy_production (xarray.DataArray): the variable that gets updated
    energy_type (string): the (new) short name of of the energy_type that will be updates
    unit (string): the units of the energy_type that should be updated.
    eur (bool): True for European sum, False for otherwise

    returns
    -------
    xarray.DataArray: input dataArray with updated attributes
    """
    short_name = "tot_" + energy_type if eur else energy_type
    what = " energy potential" if potential else " energy production"
    long_name = "European total " + energy_type + what if eur else energy_type + what

    energy_production.attrs.update(
        units=unit, short_name=short_name, long_name=long_name
    )

    return energy_production