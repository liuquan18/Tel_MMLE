# %%
# import xarray, numpy, pandas, proplot
import pandas as pd
import proplot as pplt
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import seaborn as sns
import xarray as xr
import numpy as np

# extremes
import src.extreme.extreme_ci as extreme
import src.plots.extrc_tsurf_scatter as extrc_tsurf
import src.MMLE_TEL.spatial_pattern_change as sp_change

# reimport extreme
import importlib

# warming stage
import src.warming_stage.warming_stage as warming_stage

# Stats overview
import src.plots.statistical_overview as stat_overview

# composite analysis
import src.composite.field_composite as composite

#%%
importlib.reload(extreme)
importlib.reload(sp_change)
importlib.reload(extrc_tsurf)
importlib.reload(warming_stage)
importlib.reload(stat_overview)

#%%
# class index_stats
class index_stats:
    """The class to calculate and plot the statistics of the indices"""

    def __init__(self, model, vertical_eof, fixed_pattern, standard="own") -> None:
        self.model = model
        self.vertical_eof = vertical_eof
        self.fixed_pattern = fixed_pattern
        self.standard = standard
        self.prefix = "plev_50000_" + self.fixed_pattern + "_" + self.standard

        # locations to read
        odir = "/work/mh0033/m300883/Tel_MMLE/data/" + self.model + "/"
        self.eof_result_dir = odir + "EOF_result/" + self.prefix + "_eof_result.nc"
        self.tsurf_dir = odir + "ts_processed/tsurf_mean.nc"

        # locations to save
        self.to_plot_dir = (
            "/work/mh0033/m300883/Tel_MMLE/docs/source/plots/new_standard/"
            + self.prefix
        )
        self.doc_dir = "/work/mh0033/m300883/Tel_MMLE/docs/source/"

        # read data
        self.eof_result = xr.open_dataset(self.eof_result_dir)
        self.tsurf = xr.open_dataset(self.tsurf_dir).tsurf

#%%
    # stat overview
    def stat_overview(self,levels = np.arange(-2,2.1,0.4)):
        print("ploting the statistical overview")
        first_eof = self.eof_result.eof.isel(decade = 0)
        last_eof = self.eof_result.eof.isel(decade = -1) 

        first_pc = self.eof_result.pc.isel(time = slice(0,10))
        last_pc = self.eof_result.pc.isel(time = slice(-10,None)) 

        first_fra = self.eof_result.fra.isel(decade = 0)
        last_fra = self.eof_result.fra.isel(decade = -1)


        stat_overview_fig = stat_overview.stat_overview(
            first_eof, last_eof, first_pc, last_pc, first_fra, last_fra,levels=levels
        )
        plt.savefig(self.to_plot_dir + "stat_overview.png")





    # extreme event count vs. tsurf
    def extrc_tsurf(self, plev=50000, ylim=(35, 110)):
        print("ploting the extreme event count vs. tsurf")
        tsurf_mean = self.tsurf.mean(dim = 'ens').squeeze()
        tsurf_increase = tsurf_mean - tsurf_mean[0]


        ext_counts, t_surf_mean = extrc_tsurf.decadal_extrc_tsurf(
            self.eof_result.pc, tsurf_increase
        )
        extrc_tsurf_scatter = extrc_tsurf.extCount_tsurf_scatter(
            ext_counts, t_surf_mean, ylim=ylim
        )
        plt.savefig(
            self.to_plot_dir + "extreme_count_tsurf" + f"_{(plev/100):.0f}hPa" + ".png"
        )

    
