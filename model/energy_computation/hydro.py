import xarray as xr

import config
import constants.constant as constant

def _compute_m3_capacity(capacity, height, eta=config.reservoir_efficiency):
    """ 
    Computes discharge capacities of power plants

    parameters
    ----------
    capacity (xarray.DataArray): hydro powerplant capacities
    height (xarray.DataArray): height of hydropower plants
    reservoir_efficiency (float): efficiency of hydropower plants

    returns
    -------
    m3 (xarray.DataArray): rated discharge per timestep
    """

    m3 = ((constant.MW_TO_WATT * capacity)
           /(height * constant.RHO * constant.GRAVITY * eta)
         )
    return m3

def determine_discharge_ratio(capacity, height, capacity_factor, mean_yearly_discharge):
    """ 
    computes factor-of-discharge used per gridcell based on installed capacities

    parameters
    ----------
    capacity (xarray.DataArray): hydro powerplant capacities
    height (xarray.DataArray): height of hydropower plants
    capacity_factor (float): mean yearly capcity factor of gridcell
    mean_yearly_discharge (xarray.DataArray): historic mean yearly discharge of gridcell

    returns
    -------
    discharge ratio (xarray.DataArray): ratio of discharge used in gridcell
    """
    yearly_discharge_capacity = _compute_m3_capacity(capacity, height)*(3600*24*365) 
    discharge_ratio = xr.where(
        mean_yearly_discharge>(capacity_factor*yearly_discharge_capacity),
        (capacity_factor*yearly_discharge_capacity )/mean_yearly_discharge,
        1)
    return discharge_ratio

def compute_hydro_energy_production(discharge, height, discharge_not_in_seconds = True):
    """ 
    computes factor-of-discharge used per gridcell based on installed capacities

    parameters
    ----------
    discharge(xarray.DataArray): discharge in gridcell
    height (xarray.DataArray): height of hydropower plants
    discharge_not_in_seconds (boolean): if true converted to seconds to compute Watts
    
    returns
    -------
    production (xarray.DataArray): hydropower production
    """
    production = (discharge * height * constant.RHO * constant.GRAVITY 
                  * config.reservoir_efficiency) #[W]
    if discharge_not_in_seconds:
        production = production / 3600
    production = production / constant.MW_TO_WATT #[MWh]
    return production