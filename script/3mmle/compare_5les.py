# compare the all the 5les in terms of extreme event occurence
#%%
import numpy as np
import pandas as pd
import proplot as pplt
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import seaborn as sns
import xarray as xr

# extremes
import src.extreme.extreme_ci as extreme
import src.MMLE_TEL.extrc_tsurf as extrc_tsurf
import src.MMLE_TEL.spatial_pattern_change as sp_change

# reimport extreme
import importlib

# warming stage
import src.warming_stage.warming_stage as warming_stage

# Stats overview
import src.plots.statistical_overview as stat_overview

# composite analysis
import src.composite.field_composite as composite

# %%
class compare_models:
    def __init__(self) -> None:
        self.models = [ "MK36", "CESM1_CAM5", "CanESM2", "MPI_GE"] #["GFDL_CM3",]
        self.ens_size = [30,40,50,100] # [20,30,40,50,100]
        self.odir = '/work/mh0033/m300883/Tel_MMLE/data/'
        
        self.eof_0K = [self.read_eof(model,'0K') for model in self.models]
        self.eof_4K = [self.read_eof(model,'4K') for model in self.models]

        self.extreme_0K, self.extreme_4K = self.extreme_count()
        

    def read_eof(self,model,warming_stage):
        eof_dir = self.odir + model + '/EOF_result/gph_50000_'+warming_stage+'_eof_result.nc'
        eof = xr.open_dataset(eof_dir)
        return eof

    def extreme_count(self):
        extreme_0K = [extreme.extreme_count_xr(self.eof_0K[i].pc)/self.ens_size[i] for i in range(len(self.models))]
        extreme_4K = [extreme.extreme_count_xr(self.eof_4K[i].pc)/self.ens_size[i] for i in range(len(self.models))]
        models = xr.IndexVariable('models',self.models)
        extreme_0K = xr.concat(extreme_0K,models)
        extreme_4K = xr.concat(extreme_4K,models)
        return extreme_0K,extreme_4K

    def spatial_pattern_change(self):
        pass

    def extreme_event_compare(self,mode):

        fig, axes = pplt.subplots(nrows = 1, ncols = 2)
        for i,extr_type in enumerate(['pos','neg']): # rows for pos and neg
            extreme_0K = self.extreme_0K.sel(extr_type = extr_type,mode = mode)
            extreme_4K = self.extreme_4K.sel(extr_type = extr_type,mode = mode)
            self.barplot(extreme_0K,'0K',axes[i])
            self.barplot(extreme_4K,'4K',axes[i])
        

    def barplot(self,extreme_counts,warming_stage,ax):

        true = extreme_counts.sel(confidence = 'true').values
        low = extreme_counts.sel(confidence = 'low').values
        high = extreme_counts.sel(confidence = 'high').values

        width = 0.2
        if warming_stage == '0K':
            x = 0.25 + np.arange(len(self.models))
        elif warming_stage == '4K':
            x = 0.25 + np.arange(len(self.models)) + width

        ax.errorbar(
        x = x,
        y = true,
        yerr = [(true-low),(high-true)],
        fmt = 'o',
        linewidth = 2,
        capsize = 4,
        )

        # ax.set_xticks(x + width/2, labels=self.models)
# %%
