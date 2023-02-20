# %%
# import xarray, numpy, pandas, proplot
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

# reimport extreme
import importlib

importlib.reload(extreme)

# warming stage
import src.warming_stage.warming_stage as warming_stage

# Stats overview
import src.plots.statistical_overview as stat_overview

import importlib
importlib.reload(stat_overview)

#%%
# class story line
class story_line:
    """
    generating the plots for the paper
    """
    def __init__(self,model,vertical_eof,fixed_pattern):
        self.model = model
        self.vertical_eof = vertical_eof
        self.fixed_pattern = fixed_pattern
        self.prefix = self.vertical_eof + "_" + self.fixed_pattern + "_"

        # locations
        odir = "/work/mh0033/m300883/Tel_MMLE/data/"+self.model+"/"
        self.eof_result_dir = odir + "EOF_result/"+ self.prefix + "eof_result.nc"
        self.tsurf_dir = odir + "ts_processed/tsurf_mean.nc"
        self.to_plot_dir = "/work/mh0033/m300883/Tel_MMLE/docs/source/plots/story_line/"+self.prefix

        # read eof
        self.eof_result = self.read_eof()
        self.eof = self.eof_result.eof
        self.fra = self.eof_result.fra
        self.pc = self.eof_result.pc

        # read tsurf
        self.tsurf = warming_stage.read_tsurf_fldmean(self.tsurf_dir)


        # split index into first 10 and last 10 years
        periods_pc, self.periods = warming_stage.split_period(self.pc, compare="CO2")
        self.first_pc, self.last_pc = periods_pc[0], periods_pc[1]
        self.first_eof = self.eof.sel(decade=self.periods[0]) # not time
        self.last_eof = self.eof.sel(decade=self.periods[-1])
        


        # extreme event count
        self.first_count = extreme.extreme_count_xr(self.first_pc)
        self.last_count = extreme.extreme_count_xr(self.last_pc)

    # read eof data
    def read_eof(self):
        print("reading eof result data...")
        eof_result = xr.open_dataset(self.eof_result_dir)
        return eof_result

    # statistical overview
    def stat_overview(self):
        "ploting the statistical overview"
        Fig1 = stat_overview.plot_stat_overview(self.eof_result)
        plt.savefig(self.to_plot_dir + "stat_overview.png")

    # spatial pattern change at 0K, 2K, 4K

    # extreme event count profile
    def extreme_count_profile(self):
        "ploting the extreme event count profile"
        extreme_profile = extreme.extreme_count_profile(self.first_count, self.last_count, colored=False)
        plt.savefig(self.to_plot_dir + "extreme_count_vertical_profile.png")

# *************** Figs ************************
#%%

    # plt.savefig('/work/mh0033/m300883/Tel_MMLE/docs/source/plots/story_line/Fig1.png')
# %%
# Fig 2  extreme event count profile
extreme_profile = extreme.extreme_count_profile(first_count, last_count, colored=False)
# plt.savefig('/work/mh0033/m300883/Tel_MMLE/docs/source/plots/story_line/Fig2.png')
# %%
# Fig 3  extreme event count vs. tsurf
tsurf_mean = tsurf  # .mean(dim = 'ens')
ext_counts, t_surf_mean = extrc_tsurf.decadal_extrc_tsurf(pc, tsurf)
Fig3 = extrc_tsurf.plot_extrc_tsurf(ext_counts, t_surf_mean, plev=50000)
# plt.savefig('/work/mh0033/m300883/Tel_MMLE/docs/source/plots/story_line/Fig3.png')
# %%
mpige_onepct_ind_decade = story_line("MPI_GE_onepct", "ind", "decade")
#%%
Fig1 = stat_overview.stat_overview(mpige_onepct_ind_decade.eof_result)
# %%
extreme.extreme_count_profile(mpige_onepct_ind_decade.first_count, mpige_onepct_ind_decade.last_count, colored=False,xlim = (-5,65))
# %%
