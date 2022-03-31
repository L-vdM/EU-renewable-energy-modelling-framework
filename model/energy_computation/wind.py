import xarray as xr
import numpy as np
import config
import constants.constant
import data_processing.attributes as attributes

# default values
v_cutin0= 3.5  			# cut-in wind speed [m/s]
v_rated0 = 13    		# rated wind speed [m/s]
v_cutout0 = 75#25   		# cut_out wind speed [m/s]

# =============================================================================
# Wind modules
# =============================================================================
def _wind_scale(wind_at_ref_height, hub_height, surface_roughness, reference_height):
    """
    Scales wind at reference height to wind at hub_height with a power law profile.

    parameters
    ----------
    wind_at_ref_height (xarray.DataArray): wind speed [m s**-1] 
    hub_height (float): the hub height of the wind turbines [m]
    surface_roughness (float): the roughness parameter of the surface [-]
    reference_height (float): reference height of wind speed [m].

    returns
    -------
    wind_at_hub_height (xarray.DataArray): estimated wind speed at height of the hub [m s**-1] 
    """
    wind_at_hub_height = wind_at_ref_height * (hub_height/reference_height)**surface_roughness

    return wind_at_hub_height


def _wind_power_curve(wind_at_hub_height, v_cutin, v_rated, v_cutout):
    """ 
    Computes wind energy potential with a powercurve based on wind speed at hub height.
    
    parameters
    ----------
    wind_at_hub_height (xarray.DataArray): wind speed at hub height [m s**-1]
    v_cutin (float): cut-in wind speed [m s**-1]
    v_rated (float): rated wind speed [m s**-1]
    v_cutout (float): cut-out wind speec [m s**-1]


    returns
    -------
    wind_energy_potential (xarray.DataArray): potential of wind energy at location and time [-]
    """
    wind_speed = wind_at_hub_height
    wind_energy_pot = (wind_at_hub_height**3-v_cutin**3) / (v_rated**3-v_cutin**3)
    wind_energy_pot = xr.where(wind_speed < v_cutin, 0, wind_energy_pot)    
    wind_energy_pot = xr.where(wind_speed >= v_rated, 1, wind_energy_pot)
    wind_energy_pot = xr.where(wind_speed >= v_cutout, 0, wind_energy_pot)

    return wind_energy_pot 


def compute_wind_energy_potential(
    surface_wind, 
    hub_height, 
    surface_roughness, 
    reference_height=config.height_ref, 
    v_cutin=v_cutin0, 
    v_rated=v_rated0, 
    v_cutout=v_cutout0,
    energy_type='wind',
    unit='[0-1]',
    ):
    """ 
    Computes the potential of wind energy at locations and times.
    Uses functions: wind_scale, wind_power_curve and update_energy_attributes.
    
    parameters
    ----------
    surface_wind (xarray.DataArray): wind speed at reference height [m s**-1]
    hub_height (float): the hub height of the wind turbines [m]
    surface_roughness (float): the roughness parameter of the surface [-]

    Default parameters
    ----------
    reference_height (float): reference height of wind speed [m]. Default is HEIGHT_REF
    v_cutin (float): cut-in wind speed [m s**-1]. Default is V_CUTIN
    v_rated (float): rated wind speed [m s**-1]. Default is V_RATED
    v_cutout (float): cut-out wind speec [m s**-1]. Default is V_CUTOUT

    returns
    -------
    wind_energy_pot (xarray.DataArray): potential of wind energy at location and time [-]
    """
    wind_scaled = _wind_scale(surface_wind, hub_height, surface_roughness, reference_height)
    wind_energy_pot = _wind_power_curve(wind_scaled, v_cutin, v_rated, v_cutout)
    wind_energy_pot = attributes.update_energy_attributes(
        wind_energy_pot, energy_type, unit, potential=True)
    
    return wind_energy_pot


def compute_wind_energy_production(
    wind_energy_pot,
    spatial_distribution,
    operating_time=config.time_oper,
    energy_type='wind',
    unit='MWh d**-1',
    ):
    """ 
    Computes the production of wind energy at locations and times and,
    computes the sum of wind energy production in Europe at times.
    Uses functions: _wind_scale, wind_power_curve and update_energy_attributes.
    
    parameters
    ----------
    spatial_distribution (xarray.DataArray): spatial distribution of the installed wind capacity []
    wind_energy_pot (xarray.DataArray): potential of wind energy at location and time [-]

    Default parameters
    ----------
    operatingP_time (int): operating time per day, maximum 24 [h day**-1]. Default is 'TIME_OPER'
    energy_type (string): the (new) short name of of the energy_type that will be updated. Default is 'wind'
    unit (string): the units of the energy_type that should be updated. Default is 'TWh day **-1' 

    returns
    -------
    wind_energy_prod (xarray.DataArray): production of wind energy at location and time [TWh day **-1]
    eu_wind_energy_prod (xarray.DataArray): sum of wind energy production in europe [TWh day **-1]
    """
    
    # make sure there are no rounding differences in lats between spatial and climate data
    spatial_distribution = spatial_distribution.reindex_like(
        wind_energy_pot, method='nearest', tolerance=1e-5)
    
    wind_energy_cf = wind_energy_pot * operating_time

    wind_energy_prod = wind_energy_cf * spatial_distribution

    wind_energy_prod = attributes.update_energy_attributes(
        wind_energy_prod, energy_type, unit)

    return wind_energy_prod
