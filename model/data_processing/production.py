import pandas as pd
import numpy as np 
import os
from dateutil.parser import parse

import data_processing.general as general

# def concat_country_csvs(filelist, column, timestep='1H'):
#     """
#     merge csvs after entsoe API request

#     parameters
#     ----------
#     filelist (list of strings): list of files retreived from https://transparency.entsoe.eu/
#     column (string): name of column to concat
#     timestep (string for pandas.resample): timestep to sum (get MW/hour)

#     returns
#     -------
#     pandas.DataFrame: Dataframe with MW/hour
#     """
#     df = pd.DataFrame()
#     for ifile in filelist:
#         country = os.path.basename(ifile)[0:2]
#         if country == 'EU':
#             continue   
#         dft = pd.read_csv(ifile, header=0, index_col=0, dtype='object')
#         dft = dft.loc[dft.index.notnull()]
#         cols = dft.columns
#         # convert to numeric
#         dft[cols] = dft[cols].apply(pd.to_numeric, errors='coerce')
#         # set time to utc and resampe. average MW per hour
#         dft = dft.set_index(pd.to_datetime(dft.index, utc=True)).resample(timestep).mean()
#         try:
#             series = dft[column].rename(country)
#             df = pd.concat([df, series], axis=1)
#         except:
#             continue
#     df = df.set_index(pd.to_datetime(df.index, utc=True))
#     return df

def concat_country_csvs_storage(filelist):
    """
    merge csvs of weekly data after entsoe API request. 
        row-indx=weeknr and column-indx=year

    parameters
    ----------
    filelist (list of strings): list of files retreived from https://transparency.entsoe.eu/

    returns
    -------
    pandas.DataFrame: Dataframe with week MW/hour
    """
    df = pd.DataFrame()
    for ifile in filelist:
        country = os.path.basename(ifile)[0:2] 
        dfs = pd.read_csv(ifile, index_col=0)
        # get week integers
        dfs.columns = [parse(s, fuzzy=True).year for s in dfs.columns]
        # get week integers and change to day of year
        dfs.index = [int(s)*7 for week in dfs.index for s in week.split() if s.isdigit() ]
        dfs = pd.concat([dfs[c] for c in dfs.columns], keys=dfs.columns, names=['year','week'])
        dfs = pd.DataFrame({country:dfs.apply(pd.to_numeric, errors='coerce')})
    #     list_dates = [datetime.strptime(" ".join(str(x) for x in w) + ' 0', "%Y %W %w") for w in dfs.index.to_flat_index()]
        list_dates = [general.convert_doy_to_date(w) for w in dfs.index.to_flat_index()]

        dfs.index = pd.to_datetime(list_dates)
        dfs = dfs.dropna()

        df = pd.concat([df, dfs], axis=1)
    return df


def concat_country_csvs_capacity(filelist, technology):
    """
    merge csvs of yearly installed capacities entsoe API request. 
        row-indx=weeknr and column-indx=year

    parameters
    ----------
    filelist (list of strings): list of files retreived from https://transparency.entsoe.eu/
    technology (string): tech

    returns
    -------
    pandas.DataFrame: Dataframe with week MW/hour
    """

    df = pd.DataFrame()
    for f in filelist:
        country = os.path.basename(f)[0:2]
        if country == 'EU':
            continue
        try: 
            dft = pd.read_csv(f, index_col=0)[[technology]]
            dft.index = pd.to_datetime(dft.index).tz_convert('UTC').round('1D').year
            dft = dft.rename({dft.columns[-1]: country}, axis='columns')
        except: 
            dft = pd.read_csv(f, index_col=0)
            dft.index = pd.to_datetime(dft.index).tz_convert('UTC').round('1D').year
            dft[country] = np.nan

        df = pd.concat([df,dft[[country]]], axis=1)

    df.index = pd.to_datetime(df.index,format='%Y')
    return df