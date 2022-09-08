import numpy as np
import xarray as xr
import data_processing.attributes as attributes

# =============================================================================
# Demand modules
# =============================================================================
def _weighted_area_mean(data_array, data_weights):
    """
    Computes weighted mean over a surface area longitude and latitude

    parameters
    ----------
    data_array (xarray.DataArray): array with data to be weighted wit lat and lon dimensions
    data_weights (xarray.DataArray): weights of the data with lat and lon dimensions

    returns
    -------
    weighted_average (xarray.DataArray): average over surface area (lat, lon) for each timestep
    """
    weights = data_weights / data_weights.sum()
    np.testing.assert_allclose(weights.sum().values, 1)
    weighted_average = (data_array * weights).sum(dim=["lat", "lon"], keep_attrs=True)

    return weighted_average



def LSTRmodel(temp, v):
    """
    calculates energy demand based on fit variables and (population weighted) temperature

    parameters
    ----------
    temp (xarray.DataArray): array with temperature data.
        dimemensions: (country, period, time)
    v (xarray.Set): dataset with dataarrays with variables a-f for with same dimensions as temp
        dimensions: (country, period)

    returns
    -------
    demand (xarray.DataArray)
    """
    G = 1 / (1 + np.exp(-v.e * (temp - v.f)))
    demand = (v.a + v.b * temp) * (1 - G) + (v.c + v.d * temp) * G
    return demand


def compute_demand(ds, ds_fitvalues, varin="temp", varout="demand"):
    """
    computes demand based on temperature.
    Constraints demand by maximum heating/cooling capacity based on historic entsoe demand data

    parameters
    ----------
    ds (xarray.DataSet): array with temperature data.
        dimemensions: (country, period, time)
    ds_fitvalues (xarray.DataSet): dataset with dataarrays with variables a-f and max cooling/heating
        for dimensions as temp. dimensions: (country, period)

    returns
    -------
    ds (xarray.DataSet) with demand array
    """
    ds[varout] = LSTRmodel(ds[varin], ds_fitvalues)
    # bound by highest heating and cooling capacities
    ds[varout] = xr.where(
        (ds[varout] > ds_fitvalues.heating_max) & (ds[varin] < ds_fitvalues.f),
        ds_fitvalues.heating_max,
        ds[varout],
    )
    ds[varout] = xr.where(
        (ds[varout] > ds_fitvalues.cooling_max) & (ds[varin] > ds_fitvalues.f),
        ds_fitvalues.cooling_max,
        ds[varout],
    )
    ds = attributes.set_demand_attributes(ds)
    return ds
                   
