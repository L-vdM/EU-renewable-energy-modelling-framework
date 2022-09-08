import xarray as xr
import numpy as np
import config

from constants.constant import PI, KELVIN
import data_processing.attributes as attributes

# =============================================================================
# solar modules
# =============================================================================
def _day_length(data, shift_doy=config.shift_doy):
    """
    Computes the number of daylight hours at locations and times for a given dataset

    parameters
    ----------
    dataset (xarray.DataArray(coor)): a dataset with coordinates time and lat and lon

    returns
    -------
    daylight_hours (xarray.DataArray): hours of daylight for date and location [h]
    """
    day_of_year = data.time.dt.dayofyear

    lat = data.lat

    p = 0.0
    P = np.arcsin(
        0.39795
        * np.cos(
            0.2163108
            + 2 * np.arctan(0.9671396 * np.tan(0.00860 * (day_of_year - shift_doy)))
        )
    )

    arg = (np.sin(p * PI / 180) + np.sin(lat * PI / 180) * np.sin(P)) / (
        np.cos(lat * PI / 180) * np.cos(P)
    )

    arg = xr.where(arg < -1, -1, arg)
    arg = xr.where(arg > 1, 1, arg)

    daylight_hours = 24 - 24 / PI * np.arccos(arg)

    return daylight_hours


def _solar_cell_temp(radiation, temp_at_sfc, temp_at_sfc_max, surface_wind, constants):
    """
    Computes solar cell temperature

    parameters
    ----------
    radiation (xarray.DataArray): incoming mean daily solar radiation [W m**-2]
    temp_at_sfc (xarray.DataArray): mean daily temperature at surface [degC]
    temp_at_sfc_max (xarray.DataArray): maximum daily temperature [degC]
    surface_wind (xarray.DataArray): mean daily wind at 10 m [m s**-1]
    constants (list of floats): list with 4 constants [-]

    returns
    -------
    cell_temperature (xarray.DataArray): temperature of the cell [degC]
    """
    daily_temp_2m = (temp_at_sfc + temp_at_sfc_max) / 2

    cell_temperature = (
        constants[0]
        + constants[1] * daily_temp_2m
        + constants[2] * radiation
        + constants[3] * surface_wind
    )

    return cell_temperature


def _solar_performance_ratio(cell_temperature, gamma, ref_temp):
    """
    Computes performance ratio of the

    parameters
    ----------
    cell_temperature (xarray.DataArray): temperature of the cell [degC]
    gamma (float): constant for performance ratio [-]
    ref_temp (float): reference temperature for performance ratio calculation [degC]


    returns
    -------
    performance_ratio (xarray.DataArray): performance ratio of solar cell [-]
    """
    performance_ratio = 1 + gamma * (cell_temperature - ref_temp)

    return performance_ratio


def _solar_potential(performance_ratio, radiation, standard_radiation):
    """
    Computes theoretical solar power potential

    parameters
    ----------
    performance_ratio (xarray.DataArray): performance ratio of solar cell [-]
    radiation (xarray.DataArray): incoming mean daily solar radiation [W m**-2]
    standard_radiation (int): the incoming radiation under standard test conditions [W m**-1]

    returns
    -------
    solar_energy_potential (xarray.DataArray): the daily potential for solar energy production [-]

    """
    solar_energy_potential = performance_ratio * radiation / standard_radiation

    return solar_energy_potential


def compute_solar_energy_potential(
    radiation,
    temp_at_sfc,
    temp_at_sfc_max,
    surface_wind,
    constants=config.pv_constants,
    gamma=config.gamma,
    ref_temp=config.temp_ref,
    standard_radiation=config.gstc,
    energy_type="solar",
    unit="[0-1]",
):
    """
    Computes solar energy production based on incoming radiation and performance ratio

    parameters
    ----------
    radiation (xarray.DataArray): incoming mean daily solar radiation [W m**-2]
    temp_at_sfc (xarray.DataArray): mean daily temperature at surface [Kelvin]
    temp_at_sfc_max (xarray.DataArray): maximum daily temperature [Kelvin]
    surface_wind (xarray.DataArray): mean daily wind at 10 m [m s**-1]

    predefined parameters
    ----------
    constants (list of floats): list with 4 constants for cell temperature
        calculations [-]. Default is connstants_pv
    gamma (float): constant for performance ratio [-]. Defailt is GAMMA.
    ref_temp (float): reference temperature for performance ratio
        calculation [degC]. Default is TEMP_REF
    standard_radiation (int): the incoming radiation under standard test
        conditions [W m**-1]. Default is GSTC
    energy_type (string): the (new) short name of of the energy_type that
        will be updated. Default is -solar-
    unit (string): the units of the energy_type that should be updated.
        Default is 'TWh day **-1'

    returns
    -------
    solar_energy_potential (xarray.DataArray): potential of solar energy at
        location and time [no unit]

    """
    daylight_hours = _day_length(temp_at_sfc)

    radiation_day = radiation * 24 / daylight_hours
    radiation_day = xr.where(daylight_hours == 0, 0, radiation_day)

    temp_at_sfc = temp_at_sfc - KELVIN
    temp_at_sfc_max = temp_at_sfc_max - KELVIN
    cell_temperature = _solar_cell_temp(
        radiation_day, temp_at_sfc, temp_at_sfc_max, surface_wind, constants
    )
    performance_ratio = _solar_performance_ratio(cell_temperature, gamma, ref_temp)
    solar_energy_pot = _solar_potential(
        performance_ratio, radiation_day, standard_radiation
    )

    solar_energy_pot = attributes.update_energy_attributes(
        solar_energy_pot, energy_type, unit, potential=True
    )

    return solar_energy_pot


def compute_solar_energy_production(
    solar_energy_pot, spatial_distribution, energy_type="solar", unit="MWh d**-1"
):
    """
    Computes solar energy production based on incoming radiation and performance ratio

    parameters
    ----------
    spatial_distribution (xarray.DataArray): spatial distribution of
        installed solar capacity [MW]

    predefined parameters
    ----------
    energy_type (string): the (new) short name of of the energy_type that
        will be updated. Default is 'solar'
    unit (string): the units of the energy_type that should be updated.
        Default is 'TWh day **-1'

    returns
    -------
    solar_energy_prod (xarray.DataArray): production of solar energy at
        location and time [TWh day **-1]
    eu_solar_energy_prod (xarray.DataArray): sum of solar energy production
        in europe [TWh day **-1]
    """
    # make sure there are no rounding differences in lats between spatial and climate data
    spatial_distribution = spatial_distribution.reindex_like(
        solar_energy_pot, method="nearest", tolerance=1e-5
    )

    daylight_hours = _day_length(solar_energy_pot)

    solar_energy_cf = solar_energy_pot * daylight_hours

    solar_energy_prod = solar_energy_cf * spatial_distribution
    # update attributes
    solar_energy_prod = attributes.update_energy_attributes(
        solar_energy_prod, energy_type, unit
    )

    return solar_energy_prod