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
import src.plots.extreme_plot as extreme_plot
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
importlib.reload(extreme_plot)
importlib.reload(warming_stage)
importlib.reload(stat_overview)

#%%
# class index_stats
class story_line:
    """The class to calculate and plot the statistics of the indices"""

    def __init__(
        self,
        model,
        vertical_eof,
        fixed_pattern,
        standard="temporal_ens",
        tsurf="ens_fld_year_mean",
        plev=50000,
        season="DJFM",
        tfield="DJFM",
    ) -> None:
        self.model = model
        self.vertical_eof = vertical_eof
        self.fixed_pattern = fixed_pattern
        self.standard = standard
        self.plev = plev
        self.season = season
        self.tsurf = tsurf

        self.prefix = (
            f"plev_{self.plev}_"
            + self.fixed_pattern
            + "_"
            + self.standard
            + "_"
            + self.season
        )

        # locations to read
        self.odir = "/work/mh0033/m300883/Tel_MMLE/data/" + self.model + "/"
        self.eof_result_dir = self.odir + f"EOF_result/{self.prefix}_eof_result.nc"
        self.extre_counts_500hpa_dir = (
            self.odir + f"extreme_count/{self.prefix}_extre_counts.nc"
        )
        self.extre_counts_trop_first_dir = (
            self.odir
            + f"extreme_count/troposphere_ind_decade_{self.standard}_{self.season}_first_count.nc"
        )
        self.extre_counts_trop_last_dir = (
            self.odir
            + f"extreme_count/troposphere_ind_decade_{self.standard}_{self.season}_last_count.nc"
        )
        self.tsurf_dir = self.odir + f"extreme_count/{self.tsurf}.nc"

        # locations to save
        self.to_plot_dir = (
            "/work/mh0033/m300883/Tel_MMLE/docs/source/plots/Winter_Summer/"
            + self.model
            + "_"
            + self.prefix
        )
        self.doc_dir = "/work/mh0033/m300883/Tel_MMLE/docs/source/"

        # read eof
        try:
            self.eof_result = xr.open_dataset(self.eof_result_dir)
        except FileNotFoundError:
            print("EOF result file not found.")
        try:
            self.extre_counts_500hpa = xr.open_dataset(self.extre_counts_500hpa_dir).pc
        except FileNotFoundError:
            print("500hpa extreme counts file not found.")
        try:
            self.extre_counts_trop_first = xr.open_dataset(
                self.extre_counts_trop_first_dir
            ).pc
        except FileNotFoundError:
            print("Troposphere first extreme counts file not found.")
        try:
            self.extre_counts_trop_last = xr.open_dataset(
                self.extre_counts_trop_last_dir
            ).pc
        except FileNotFoundError:
            print("Troposphere last extreme counts file not found.")
        try:
            self.tsurf = self.read_tsurf()
        except FileNotFoundError:
            print("tsurf file not found.")

    #%%
    # read tsurf
    def read_tsurf(self):
        tsurf = xr.open_dataset(self.tsurf_dir)

        # read the array either tusrf, ts, or tas
        try:
            tsurf_arr = tsurf.tsurf.squeeze()
        except AttributeError:
            tsurf_arr = tsurf.ts.squeeze()
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
        # extreme counts of first and last 10 decades

        print("ploting the extreme event count profile")
        extreme_profile = extreme_plot.extreme_count_profile(
            self.extre_counts_trop_first,
            self.extre_counts_trop_last,
            colored=False,
            **kwargs,
        )
        plt.savefig(
            self.to_plot_dir.replace(self.prefix, "")
            + self.season
            + "_"
            + "extreme_count_vertical_profile.png"
        )
        # slightly different path for the fig.

    # extreme event count vs. tsurf
    def extrc_tsurf(self, ylim=(35, 110)):

        extrc_tsurf_scatter = extreme_plot.extCount_tsurf_scatter(
            self.extre_counts_500hpa, self.tsurf, ylim=ylim
        )
        plt.savefig(self.to_plot_dir + f"_extreme_count_tsurf.png", dpi=300)

    # composite analysis of surface temperature in terms of different extreme events
    def composite_analysis(self, reduction="mean", tfield = 'same',level_bound = None, levels_NAO = np.arange(-3,3.1,0.5), levels_EA = np.arange(-2,2.1,0.5)):
        """
        tfield can be 'same' or 'next'
        """
        if tfield == 'same':
            tsurf_season = self.season
        else:
            if tfield == 'next':
                tsurf_season = self.get_next_season()
            elif tfield == 'DJF' or tfield == 'JJA' or tfield == 'SON' or tfield == 'MAM' or tfield == 'JJAS' or tfield == 'DJFM':
                tsurf_season = tfield
            else:
                print("wrong input of tfield")
                tsurf_season = None
        first_com_tsurf_dir = f"{self.odir}composite/{self.prefix}_{tsurf_season}_first_composite.nc"
        last_com_tsurf_dir = f"{self.odir}composite/{self.prefix}_{tsurf_season}_last_composite.nc"

        print("ploting the composite analysis of surface temperature")
        print(f" reading the composite data of {self.prefix}_{tsurf_season}")

        first_var = xr.open_dataset(first_com_tsurf_dir)
        last_var = xr.open_dataset(last_com_tsurf_dir)

        try:
            first_var = first_var.tsurf.squeeze()
        except AttributeError:
            try:
                first_var = first_var.ts.squeeze()
            except AttributeError:
                first_var = first_var.tas.squeeze()

        try:
            last_var = last_var.tsurf.squeeze()
        except AttributeError:
            try:
                last_var = last_var.ts.squeeze()
            except AttributeError:
                last_var = last_var.tas.squeeze()

        temp_NAO = composite.composite_plot(first_var, last_var, "NAO",level_bound = level_bound,levels=levels_NAO)
        plt.savefig(
            self.to_plot_dir + f"_{self.prefix}_{tsurf_season}_composite_tsurf_NAO.png",
            dpi=300,
        )

        # temp_EA = composite.composite_plot(first_var, last_var, "EA",level_bound = level_bound,levels=levels_EA)
        # plt.savefig(
        #     self.to_plot_dir + f"_{season_field_dir}_composite_tsurf_EA.png",
        #     dpi=300,
        # )

    def get_next_season(self):
        if len(self.season) == 4: # DJFM and JJAS
            seasons = ['DJFM','MAM','JJAS','SON']
        elif len(self.season) == 3: # DJF and JJA
            seasons = ['DJF','MAM','JJA','SON']

        index = seasons.index(self.season)
        next_index = (index + 1)% len(seasons)
        next_season = seasons[next_index]
        return next_season

    # plot all
    def plot_all(self):
        self.stat_overview()
        self.NAO_EA_hist2d()
        self.extreme_count_profile()
        self.extrc_tsurf()

    # write the above plots into a markdown file
    def write_doc(self):
        """create the md file for the plots"""
        relative_plot_dir = (
            "plots/Winter_Summer/"
            + self.model
            + "_plev_50000_"
            + self.fixed_pattern
            + "_"
            + self.standard
            + "_"
            + self.season
            + "_"
        )
        print("creating the markdown file for the plots")

        next_season = self.get_next_season()

        with open(
            self.doc_dir
            + self.model
            + "_"
            + self.fixed_pattern
            + "_"
            + self.season
            + "_index_stats.md",
            "w",
        ) as f:
            f.write(f"# Statistics of the indices in {self.season}\n")

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

            f.write("## extreme event count profile\n")
            f.write(
                f"![extreme event count profile](plots/Winter_Summer/{self.model}_{self.season}_extreme_count_vertical_profile.png)\n"
            )

            f.write("## extreme event count  vs. tsurf\n")
            f.write(
                f"![extreme event count vs. tsurf]({relative_plot_dir}extreme_count_tsurf.png)\n"
            )

            f.write(
                "## composite analysis of surface temperature in terms of different extreme events\n"
            )
            f.write(
                f"![composite analysis of surface temperature in terms of different extreme events]({relative_plot_dir}{self.season}_composite_tsurf_NAO.png)\n"
            )
            f.write(
                f"![composite analysis of surface temperature in terms of different extreme events]({relative_plot_dir}{self.season}_composite_tsurf_EA.png)\n"
            )

            f.write(
                f"![composite analysis of surface temperature in next season]({relative_plot_dir}{next_season}_composite_tsurf_NAO.png)\n"
            )
            f.write(
                f"![composite analysis of surface temperature in next season]({relative_plot_dir}{next_season}_composite_tsurf_EA.png)\n"
            )

