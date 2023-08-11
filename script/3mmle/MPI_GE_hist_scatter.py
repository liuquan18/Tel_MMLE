#%%
# imports
import xarray as xr
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import proplot as pplt
import seaborn as sns
from pandas.tseries.offsets import DateOffset
import warnings

# functions to plot
import src.plots.vertical_profile as profile_plots
import src.plots.PDF as pdf_plots
import src.plots.plot_violin as violin_plots
import src.plots.spatial_distribution_plot as spatial_dis_plots
import src.plots.return_period as RP_plots
import src.plots.composite_spatial_pattern as composite_spatial_pattern
import src.plots.composite_plot as composite_plot
import src.warming_stage.warming_stage as warming_stage

import src.extreme.period_pattern_extreme as extreme
import src.EVT.return_period as EVT
import src.composite.field_composite as composite
import src.html.create_md as create_md
import src.Teleconnection.tools as tools
import src.MMLE_TEL.spatial_pattern_change as sp_change
import src.MMLE_TEL.extrc_tsurf as extrc_tsurf
import warnings

warnings.filterwarnings("ignore")

#%%
import importlib
importlib.reload(extrc_tsurf)

# %%
# read eof
eof_result = xr.open_dataset("/work/mh0033/m300883/Tel_MMLE/data/MPI_GE/EOF_result/gph_50000_decade_eof_result.nc")

#%%
# read temperature
ts_mean = xr.open_dataset("/work/mh0033/m300883/Tel_MMLE/data/MPI_GE/ts_processed/ens_fld_year_mean.nc")
ts_mean = ts_mean.mean(dim = 'ens')
ts_mean = ts_mean - ts_mean.isel(time = slice(0,10)).mean(dim = 'time')

# %%
extr_count, ts_mean = extrc_tsurf.decadal_extrc_tsurf(
            eof_result.pc, ts_mean
        )
# %%
# assign new coordinate called 'plev', with value of 50000 to ts_mean
ts_mean = ts_mean.assign_coords(plev = 50000)
# same to extr_count
extr_count = extr_count.assign_coords(plev = 50000)
#%%
fig = extrc_tsurf.extCount_tsurf_scatter(extr_count, ts_mean,plev = 50000)
# %%
fig, axes = pplt.subplots(nrows=2, ncols=2, figwidth=8, span=False, share=False,ylim = (470,530))



# %%
