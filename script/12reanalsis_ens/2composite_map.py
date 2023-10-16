# %%
import xarray as xr
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib as mpl
import proplot as pplt
import seaborn as sns
import cartopy.crs as ccrs
import matplotlib.ticker as ticker
from matplotlib import lines as mlines
import matplotlib.ticker as ticker



from matplotlib.ticker import MultipleLocator, FormatStrFormatter, MaxNLocator
import matplotlib.patches as mpatches


import src.plots.composite_plot as composite_plot
import src.plots.extreme_plot as extplt
import src.plots.statistical_overview as stat_overview
import src.plots.utils as utils
import src.obs.era5_extreme_change as era5_extreme_change
#%%
import importlib
importlib.reload(composite_plot)
# %%
first_ts = xr.open_dataset("/work/mh0033/m300883/Tel_MMLE/data/ERA5_allens/composite/first_composite.nc")
last_ts = xr.open_dataset("/work/mh0033/m300883/Tel_MMLE/data/ERA5_allens/composite/last_composite.nc")
# %%
first_ts = first_ts.__xarray_dataarray_variable__
last_ts = last_ts.__xarray_dataarray_variable__
# %%
composite_plot.composite_plot(first_ts, last_ts, 'NAO', levels=np.arange(-1.5, 1.6, 0.3))
plt.savefig("/work/mh0033/m300883/Tel_MMLE/docs/source/plots/supplyment/composite_ERA5_allens.png", dpi=300, bbox_inches='tight')
# %%
