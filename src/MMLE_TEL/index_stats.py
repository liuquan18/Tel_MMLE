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
import src.plots.NAO_EA_hist2d as hist2d


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

    def __init__(
        self, model, vertical_eof, fixed_pattern, standard="temporal_ens", tsurf = 'ens_fld_year_mean', plev = 50000,season = 'DJFM'
    ) -> None:
        self.model = model
        self.vertical_eof = vertical_eof
        self.fixed_pattern = fixed_pattern
        self.standard = standard
        self.plev = plev
        self.season = season
        self.prefix = f"plev_{self.plev}_" + self.fixed_pattern + "_" + self.standard + "_" + self.season


        # locations to read
        self.odir = "/work/mh0033/m300883/Tel_MMLE/data/" + self.model + "/"
        self.eof_result_dir = self.odir + "EOF_result/" + self.prefix + "_eof_result.nc"
        self.tsurf_dir = self.odir + "ts_processed/" + tsurf + ".nc"
        self.field_tsurf_dir = self.odir + "ts/" + tsurf + ".nc"

        # locations to save
        self.to_plot_dir = (
            "/work/mh0033/m300883/Tel_MMLE/docs/source/plots/new_standard/"
            + self.model
            + "_"
            + self.prefix
        )
        self.doc_dir = "/work/mh0033/m300883/Tel_MMLE/docs/source/"

        # read eof
        self.eof_result = xr.open_dataset(self.eof_result_dir)
        self.tsurf = self.read_tsurf()

    #%%
    # read tsurf
    def read_tsurf(self):
        tsurf = xr.open_dataset(self.tsurf_dir)

        # read the array either tusrf, ts, or tas
        try:
            tsurf_arr = tsurf.tsurf.squeeze()
        except AttributeError:
            try:
                tsurf_arr = tsurf.ts.squeeze()
            except AttributeError:
                tsurf_arr = tsurf.tas.squeeze()
        # change the temp time into datetime64
        try:
            tsurf_arr["time"] = tsurf_arr.indexes["time"].to_datetimeindex()
        except AttributeError:
            pass
        return tsurf_arr

    #%%
    # stat overview
    def stat_overview(self, levels=np.arange(-2, 2.1, 0.4)):
        print("ploting the statistical overview")
        first_pc = self.eof_result.pc.isel(time=slice(0, 10))
        last_pc = self.eof_result.pc.isel(time=slice(-10, None))
        first_eof = None
        last_eof = None
        first_fra = None
        last_fra = None

        if self.fixed_pattern == "decade":
            first_eof = self.eof_result.eof.isel(decade=0)
            last_eof = self.eof_result.eof.isel(decade=-1)

            first_fra = self.eof_result.fra.isel(decade=0)
            last_fra = self.eof_result.fra.isel(decade=-1)
        elif self.fixed_pattern == "all":
            first_eof = self.eof_result.eof
            last_eof = None

            first_fra = self.eof_result.fra
            last_fra = None



        stat_overview_fig = stat_overview.stat_overview(
            first_pc, last_pc, first_fra, last_fra, first_eof, last_eof, levels=levels
        )
        plt.savefig(self.to_plot_dir + "_stat_overview.png", dpi=300)

    # 2d hist of NAO and EA in the first10 and last10 decades
    def NAO_EA_hist2d(self, bins=20, levels=np.arange(0, 0.31, 0.03)):
        print("ploting the 2d hist of NAO and EA in the first10 and last10 decades")
        hist2d.NAO_EA_hist2d(self.eof_result.pc, bins=bins, levels=levels)
        plt.savefig(self.to_plot_dir + "_NAO_EA_hist2d.png", dpi=300)

    # plot the vertical profile of the extreme counts
    def extreme_count_profile(self, **kwargs):
        print("reading the eof result with whole troposphere")
        ci = kwargs.pop("ci", "AR1")
        # read the eof result with whole troposphere
        eof_result_trop = xr.open_dataset(
            "/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_onepct/EOF_result/"
            + "troposphere_"
            + self.vertical_eof
            + "_"
            + self.fixed_pattern
            + "_"
            + self.standard
            + "_eof_result.nc"
        )

        # pcs of first and last 10 decades
        first_pc = eof_result_trop.pc.isel(time=slice(0, 10))
        last_pc = eof_result_trop.pc.isel(time=slice(-10, None))

        print("calculating the extreme event count")
        # extreme counts of first and last 10 decades
        first_count = extreme.extreme_count_xr(first_pc, ci=ci)
        last_count = extreme.extreme_count_xr(last_pc, ci=ci)

        print("ploting the extreme event count profile")
        extreme_profile = extreme.extreme_count_profile(
            first_count, last_count, colored=False, **kwargs
        )
        plt.savefig(
            self.to_plot_dir[:-44]
            + f"{self.model}_{ci}_extreme_count_vertical_profile.png"
        )
        # slightly different path for the fig.

    # extreme event count vs. tsurf
    def extrc_tsurf(self, ylim=(35, 110), ci="AR1"):

        # check if the file exists
        try:
            ds = xr.open_dataset(
                self.odir + 'extreme_count/extre_counts_tsurf.nc'
            )
            ext_counts = ds.extreme_counts
            t_surf_mean = ds.tsurf

        except FileNotFoundError:
            print("ploting the extreme event count vs. tsurf")
            try:
                tsurf_mean = self.tsurf.mean(dim="ens").squeeze()
            except ValueError:
                tsurf_mean = self.tsurf
            tsurf_increase = tsurf_mean - tsurf_mean[0]

            ext_counts, t_surf_mean = extrc_tsurf.decadal_extrc_tsurf(
                self.eof_result.pc, temp = tsurf_increase, ci=ci
            )
            
        extrc_tsurf_scatter = extrc_tsurf.extCount_tsurf_scatter(
            ext_counts, t_surf_mean, ylim=ylim
        )
        plt.savefig(self.to_plot_dir + f"_{ci}_extreme_count_tsurf.png", dpi=300)

    # composite analysis of surface temperature in terms of different extreme events
    def composite_analysis(self, reduction="mean", **kwargs):
        print("ploting the composite analysis of surface temperature")
        print(" reading the tsurf data...")
        var_data = xr.open_dataset(self.field_tsurf_dir + "all_ens_tsurf.nc").tsurf
        var_data = var_data - var_data.mean(dim="ens")

        first_index = self.eof_result.pc.isel(time=slice(0, 10))
        last_index = self.eof_result.pc.isel(time=slice(-10, None))

        print(" compositing the tsurf data...")
        first_var = composite.Tel_field_composite(
            first_index, var_data, threshold=1.5, reduction=reduction, **kwargs
        )
        last_var = composite.Tel_field_composite(
            last_index, var_data, threshold=1.5, reduction=reduction, **kwargs
        )

        temp_NAO = composite.composite_plot(first_var, last_var, "NAO")
        plt.savefig(
            self.to_plot_dir + f"_composite_tsurf_NAO.png",
            dpi=300,
        )

        temp_EA = composite.composite_plot(first_var, last_var, "EA")
        plt.savefig(
            self.to_plot_dir + f"_composite_tsurf_EA.png",
            dpi=300,
        )

    # plot all
    def plot_all(self):
        self.stat_overview()
        self.NAO_EA_hist2d()
        self.extrc_tsurf()

    # write the above plots into a markdown file
    def write_doc(self):
        """create the md file for the plots"""
        relative_plot_dir = (
            "plots/new_standard/"
            + self.model
            + "_plev_50000_"
            + self.fixed_pattern
            + "_"
            + self.standard
            + "_"
        )
        print("creating the markdown file for the plots")
        with open(
            self.doc_dir + self.model + "_" + self.fixed_pattern + "_index_stats.md",
            "w",
        ) as f:
            f.write("# Statistics of the indices\n")

            f.write("This file contains the statistics of the indices\n")

            f.write("## statistical overview\n")
            f.write(
                "The first EOF and PC of the first 10 decades and the last 10 decades are shown below\n"
            )
            f.write(f"![statistical overview]({relative_plot_dir}stat_overview.png)\n")

            f.write("## 2d hist of NAO and EA in the first10 and last10 decades\n")
            f.write(
                f"![2d hist of NAO and EA in the first10 and last10 decades]({relative_plot_dir}NAO_EA_hist2d.png)\n"
            )

            f.write("## extreme event count profile with AR(1)\n")
            f.write(
                f"![extreme event count profile](plots/new_standard/{self.model}_AR1_extreme_count_vertical_profile.png)\n"
            )
            f.write("## extreme event count profile with bootstrap\n")
            f.write(
                f"![extreme event count profile](plots/new_standard/{self.model}_bootstrap_extreme_count_vertical_profile.png)\n"
            )

            f.write("## extreme event count AR(1) vs. tsurf\n")
            f.write(
                f"![extreme event count vs. tsurf]({relative_plot_dir}AR1_extreme_count_tsurf.png)\n"
            )

            f.write("## extreme event count bootstrap vs. tsurf\n")
            f.write(
                f"![extreme event count vs. tsurf]({relative_plot_dir}bootstrap_extreme_count_tsurf.png)\n"
            )

            f.write(
                "## composite analysis of surface temperature in terms of different extreme events\n"
            )
            f.write(
                f"![composite analysis of surface temperature in terms of different extreme events]({relative_plot_dir}composite_tsurf_NAO.png)\n"
            )
            f.write(
                f"![composite analysis of surface temperature in terms of different extreme events]({relative_plot_dir}composite_tsurf_EA.png)\n"
            )
