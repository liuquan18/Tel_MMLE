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

    slope = slope 
    confidence_interval = confidence_interval 
    
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

models_all = ["MPI_GE", "CanESM2", "CESM1_CAM5", "MK36", "GFDL_CM3"]
colors_model = {"MPI_GE": "C1", "CanESM2": "tab:purple", "CESM1_CAM5": "tab:blue", "MK36": "tab:green", "GFDL_CM3": "yellow"}

models_legend = {"MPI_GE": "MPI-GE (100)", "CanESM2": "CanESM2 (50)", "CESM1_CAM5": "CESM1-CAM5 (40)", "MK36": "MK3.6 (30)", "GFDL_CM3": "GFDL-CM3 (20)"}

fig, ax = plt.subplots()
for model in models_all:
    data = flow_regimes_df[flow_regimes_df.model == model]
    color = colors_model[model]
    label = models_legend[model]

    sns.scatterplot(data=data, x='jet_loc', y='GB', ax=ax, color=color, label=label, s=100)
    ax.plot([data.jet_loc_CI.iloc[0][0], data.jet_loc_CI.iloc[0][1]], [data.GB.iloc[0], data.GB.iloc[0]], color=color)
    ax.plot([data.jet_loc.iloc[0], data.jet_loc.iloc[0]], [data.GB_CI.iloc[0][0], data.GB_CI.iloc[0][1]], color=color)



ax.legend(loc = 'lower right', frameon=False)
plt.show()

# %%
