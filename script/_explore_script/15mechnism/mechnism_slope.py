#%%
import xarray as xr
import pandas as pd
import seaborn as sns
import glob
import matplotlib.pyplot as plt
import os
import numpy as np
from scipy import stats
# %%
# read jet_loc and GB clim
def read_clim(model = 'MPI_GE', extreme = False):
    dir = f"/work/mh0033/m300883/Tel_MMLE/data/{model}/mechnisms/"

    if extreme:
        jet_loc_clim = xr.open_dataset(dir + 'jet_north_decade.nc').jet_loc
        GB_clim = xr.open_dataset(dir + 'GB_pos_decade.nc').GB
        return jet_loc_clim, GB_clim
    
    jet_loc_clim = xr.open_dataset(dir + 'jet_loc_clim.nc').jet_loc
    GB_clim = xr.open_dataset(dir + 'GB_clim.nc').GB
    return jet_loc_clim, GB_clim


# %%
def calculate_slope_and_confidence_interval(var_clim):

    # standardize the data
    # var_clim = (var_clim - var_clim.mean()) / var_clim.std()


    # Calculate the slope and intercept using linear regression
    x = np.arange(len(var_clim.time))
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, var_clim)
    
    # Calculate the 5%-95% confidence interval
    confidence_interval = 1.96 * std_err

    # into decade (*10)
    slope = slope * 10
    confidence_interval = confidence_interval * 10
    
    return slope, (slope - confidence_interval, slope + confidence_interval)

# %%
jet_slopes = {}
jet_confidence_intervals = {}

GB_slopes = {}
GB_confidence_intervals = {}
# %%
for model in ["MPI_GE", "CanESM2", "CESM1_CAM5", "GFDL_CM3", "MK36"]:
    jet_loc_clim, GB_clim = read_clim(model, False)
    jet_slope, jet_confidence_interval = calculate_slope_and_confidence_interval(jet_loc_clim)
    jet_slopes[model] = jet_slope
    jet_confidence_intervals[model] = jet_confidence_interval

    GB_slope, GB_confidence_interval = calculate_slope_and_confidence_interval(GB_clim)
    GB_slopes[model] = GB_slope
    GB_confidence_intervals[model] = GB_confidence_interval
# %%
# into dataframes
flow_regimes_df = pd.DataFrame(jet_slopes, index = ['jet_loc']).T 
# %%
flow_regimes_df['jet_loc_CI'] = flow_regimes_df.index.map(jet_confidence_intervals)
# %%
flow_regimes_df['GB'] = flow_regimes_df.index.map(GB_slopes)
# %%
flow_regimes_df['GB_CI'] = flow_regimes_df.index.map(GB_confidence_intervals)
#%%
flow_regimes_df = flow_regimes_df.reset_index().rename(columns = {'index':'model'})
# %%
fig, ax = plt.subplots()
# plot scatter of x-jet_loc and y-GB
sns.scatterplot(data = flow_regimes_df, x = 'jet_loc', y = 'GB', ax = ax, hue = 'model')
# plot the confidence interval
for i, row in flow_regimes_df.iterrows():
    ax.plot([row.jet_loc_CI[0], row.jet_loc_CI[1]], [row.GB, row.GB], color = 'gray')
    ax.plot([row.jet_loc, row.jet_loc], [row.GB_CI[0], row.GB_CI[1]], color = 'gray')
# %%
