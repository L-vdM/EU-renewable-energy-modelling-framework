# import pandas as pd

# def _add_netto_demand(df):
#     """
#     add netto demand to dataset

#     parameters
#     ----------
#     df (pandas.DataFrame): dataset with demand and wind and solar production and run-f-ruver

#     returns
#     -------
#     mdf (pandas.DataFrame): with netto demand
#     """
#     df['netto_demand'] = (df['demand']-df['wind_solar'])
#     df.loc[df['netto_demand']<0, 'netto_demand'] = 0
#     return df

# def read_xarrays_to_df(demandf, inflowf, wind_and_solar, mapping):
#     """
#     read netcdf demand, inflow and wind and solar production and join in pandas dataframe

#     parameters
#     ----------
#     demandf,inflowf,wind_and_solar (string): file names of netcdf with productions and demand


#     returns
#     -------
#     df (pandas.DataSet): with demand, production and netto demand
#     """
#     ### open demand
#     demand = xr.open_dataset(demandf)
#     # demand has double timestamps from weekend/weekday data
#     demand = demand.dropna(dim='time')
#     # dropweighetd temperature data from demand
#     demand = demand.drop(['temp'], axis=1)

#     ### open inflow [MWh]
#     inflow = xr.open_dataset(inflowf)
#     ### open productio.
#     prod = xr.open_dataset(wind_and_solar)
#     # sum all vars in production (wind and solar)
#     vars_to_sum = list(prod.keys())
#     prod_total = prod[vars_to_sum].to_array().sum("variable").to_dataset(name='wind_solar')

#     ### concat data
#     dc = demand.drop(['period']).to_dataframe()
#     ic = inflow.to_dataframe()
#     pc = prod_total.to_dataframe()
#     df = pd.concat([dc, ic, pc], axis=1)
#     # change all units the MWh
#     df['demand'] = df['demand']*1e3

#     return _add_netto_demand(df)




