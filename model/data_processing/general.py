from datetime import datetime, timedelta
import xarray as xr
import pandas as pd
import numpy as np

# from sklearn.metrics import mean_squared_error
# from scipy import stats

# def wavg(group, avg_name, weight_name):
#     """ http://stackoverflow.com/questions/10951341/pandas-dataframe-aggregate-function-using-multiple-columns
#     In rare instance, we may not have weights, so just return the mean. Customize this if your business case
#     should return otherwise.
#     """
#     d = group[avg_name]
#     w = group[weight_name]
#     try:
#         return (d * w).sum() / w.sum()
#     except ZeroDivisionError:
#         return d.mean()


def convert_doy_to_date(doy):
    """
    converts year and day of year to datetime object

    parameters
    ----------
    doy (tuple or list): (year, dayofyear)

    returns
    -------
    datetime.datetime: datetime object with year,month,day
    """
    return datetime(doy[0], 1, 1) + timedelta(doy[1] - 1)


# def xr_lat_lon_to_flat_df(ds, var):
#     df = ds[var].to_dataframe()
#     df.reset_index(inplace=True, level = ['lat', 'lon'])
#     df = df.dropna()
#     df['loc'] = df[['lat', 'lon']].apply(tuple, axis=1)

#     try:
#         df = df.drop(['lat', 'lon', 'country'], axis=1)
#     except:
#         df = df.drop(['lat', 'lon'], axis=1)
#     df = df.reset_index().set_index(['loc',"time"]).unstack(level=0)
#     return df

# def is_outlier(points, thresh=3.5):
#     """
#     Returns a boolean array with True if points are outliers and False
#     otherwise.

#     Parameters:
#     -----------
#         points : An numobservations by numdimensions array of observations
#         thresh : The modified z-score to use as a threshold. Observations with
#             a modified z-score (based on the median absolute deviation) greater
#             than this value will be classified as outliers.

#     Returns:
#     --------
#         mask : A numobservations-length boolean array.

#     References:
#     ----------
#         Boris Iglewicz and David Hoaglin (1993), "Volume 16: How to Detect and
#         Handle Outliers", The ASQC Basic References in Quality Control:
#         Statistical Techniques, Edward F. Mykytka, Ph.D., Editor.
#     """
#     if len(points.shape) == 1:
#         points = points[:,None]
#     median = np.median(points, axis=0)
#     diff = np.sum((points - median)**2, axis=-1)
#     diff = np.sqrt(diff)
#     med_abs_deviation = np.median(diff)

#     modified_z_score = 0.6745 * diff / med_abs_deviation

#     return modified_z_score > thresh

# def get_flat_pandaframe_of_2D_dataSet(ds, var):
#     df = ds[var].to_dataframe().reset_index()
#     df.replace('', np.nan, inplace=True)
#     df = df.dropna()
#     df = df.set_index(var)


# def compute_weighted_var(ds, da_weights):
#     ds = (ds*da_weights)
#     ds = ds.sum(dim=['lat', 'lon'])
#     return ds

# def find_linear_intersect(a, b, c, d):
#     x = (a-c)/ (d-b)
#     y = a + b*x
#     return x, y

# def get_linear_part(S, L, xvar='temp', yvar='fdemand'):
#     E_best = np.inf
#     k_best = 0
#     q_best = 0
#     for k in S.index:
#         for q in range(k+2, S.index[-1]):
#                 if (S[xvar][q]-S[xvar][k])>L and (S[xvar][q-1]-S[xvar][q]<L or S[xvar][q-1]+S[xvar][q]<L):
#                     slope, intercept, r_value, p_value, std_err = stats.linregress(
#                         S[k:q][xvar], S[k:q][yvar])

#                     y_fit = intercept + S[k:q][xvar] * slope
#                     rms = mean_squared_error(S[k:q][yvar], y_fit, squared=False)
#                     E = rms
#                     if E<E_best:
#                         E_best = E
#                         k_best = k
#                         q_best = q
#     return k_best,q_best

# def concat_df_over_variable(df, varname, dimname):
#     var_values = np.unique(df[varname].values)
#     var_values = var_values[np.where(var_values!='')]
#     dfc = xr.concat([df.where(df[varname]==v) for v in var_values], dimname)
#     dfc[dimname] = var_values
#     return dfc


def sum_per_region(region_mask, ds):
    dfc = (region_mask.mask * ds).sum(dim=["lat", "lon"])
    return dfc