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
# class story line
class story_line:
    """
    generating the plots for the paper
    """

    def __init__(self, model, vertical_eof, fixed_pattern):
        self.model = model
        self.vertical_eof = vertical_eof
        self.fixed_pattern = fixed_pattern
        self.prefix = self.vertical_eof + "_" + self.fixed_pattern + "_"

        # locations
        odir = "/work/mh0033/m300883/Tel_MMLE/data/" + self.model + "/"
        self.eof_result_dir = odir + "EOF_result/" + self.prefix + "own_eof_result.nc"
        self.tsurf_dir = odir + "ts_processed/tsurf_mean.nc"
        self.ts_fullfield_dir = (
            odir + "ts/"
        )  # the location for the full field temperature data
        self.atg_dir = odir + "ts_processed/tsurf_anom_gradient.nc"

        ## locations to save
        self.to_plot_dir = (
            "/work/mh0033/m300883/Tel_MMLE/docs/source/plots/story_line/" + self.prefix
        )
        self.doc_dir = "/work/mh0033/m300883/Tel_MMLE/docs/source/"

        # read eof
        self.eof_result = self.read_eof()
        self.eof = self.eof_result.eof
        self.fra = self.eof_result.fra
        self.pc = self.eof_result.pc

        # read tsurf
        print("reading tsurf data...")
        self.tsurf = warming_stage.read_tsurf_fldmean(self.tsurf_dir)
        print("reading tsurf gradient data...")
        self.atg = warming_stage.read_tsurf_fldmean(self.atg_dir)
        self.warming_periods = warming_stage.temp_period(self.tsurf)

        # split index into first 10 and last 10 years
        periods_pc, self.periods = warming_stage.split_period(self.pc, compare="CO2")
        self.first_pc, self.last_pc = periods_pc[0], periods_pc[1]
        self.first_eof = self.eof.sel(decade=self.periods[0])  # not time
        self.last_eof = self.eof.sel(decade=self.periods[-1])

        # extreme event count
        self.first_count = extreme.extreme_count_xr(self.first_pc)
        self.last_count = extreme.extreme_count_xr(self.last_pc)

    # read eof data
    def read_eof(self):
        print("reading eof result data...")
        eof_result = xr.open_dataset(self.eof_result_dir)
        return eof_result

    # read temperature data to do the composite analysis
    def read_var(self):
        """
        read the tsurf or wind, precipitation for composite analysis
        data are stored in 3rdPanel/data/
        """
        print("reading the var to do composite...")
        # MPI is different format...
        if self.model == "MPI_GE_onepct":
            var_data = xr.open_dataset(self.ts_fullfield_dir + "all_ens_tsurf.nc").tsurf
        elif self.model == "MPI_GE":
            var_data = xr.open_mfdataset(
                self.ts_fullfield_dir + "*.nc", combine="nested", concat_dim="ens"
            ).tsurf
        else:
            var_data = xr.open_mfdataset(
                self.ts_fullfield_dir + "*.nc", combine="nested", concat_dim="ens"
            )
            var_data = var_data.ts
            var_data = var_data.rename({"plev": "plev"})
            var_data["time"] = var_data.indexes["time"].to_datetimeindex()
        var_data = var_data - var_data.mean(dim="ens")  # remove the ensemble mean

        return var_data

    # statistical overview
    def stat_overview(self):
        print("ploting the statistical overview")
        Fig1 = stat_overview.stat_overview(self.eof_result)
        plt.savefig(self.to_plot_dir + "stat_overview.png")

    # spatial pattern change at 0K, 2K, 4K
    def spatial_pattern_change(self, plev=50000, **kwargs):
        print("ploting the spatial pattern change at 0K, 2K, 4K")
        Fig2 = sp_change.spatial_pattern_change_decade(
            self.warming_periods, self.eof, self.fra, **kwargs
        )
        plt.savefig(
            self.to_plot_dir
            + "spatial_pattern_change"
            + f"_{(plev/100):.0f}hPa"
            + ".png"
        )

    # extreme event count profile
    def extreme_count_profile(self, **kwargs):
        print("ploting the extreme event count profile")
        extreme_profile = extreme.extreme_count_profile(
            self.first_count, self.last_count, colored=False, **kwargs
        )
        plt.savefig(self.to_plot_dir + "extreme_count_vertical_profile.png")

    # extreme event count vs. tsurf
    def extrc_tsurf(self, plev=50000, ylim=(5, 65)):
        print("ploting the extreme event count vs. tsurf")
        tsurf_mean = self.tsurf
        tsurf_increase = tsurf_mean - tsurf_mean[0]
        ext_counts, t_surf_mean = extrc_tsurf.decadal_extrc_tsurf(
            self.pc, tsurf_increase
        )
        Fig3 = extrc_tsurf.extCount_tsurf_scatter(
            ext_counts, t_surf_mean, plev=plev, ylim=ylim
        )
        plt.savefig(
            self.to_plot_dir + "extreme_count_tsurf" + f"_{(plev/100):.0f}hPa" + ".png"
        )

    # extreme event count vs. arctic-tropical gradient
    def extrc_atg(self, plev=50000, ylim=(-5, 65)):
        print("ploting the extreme event count vs. arctic-tropical gradient")
        atg = self.atg
        ext_counts, atg_dec = extrc_tsurf.decadal_extrc_tsurf(self.pc, atg)
        extc_atg_scatter = extrc_tsurf.extCount_tsurf_scatter(
            ext_counts,
            atg_dec,
            plev=plev,
            ylim=ylim,
            xlim=(-7, 7),
            xlabel="tropic-arctic gradient",
        )
        plt.savefig(
            self.to_plot_dir + "extreme_count_atg" + f"_{(plev/100):.0f}hPa" + ".png"
        )

    # composite analysis of surface temperature in terms of different extreme events
    def composite_analysis(self, plev=50000, **kwargs):
        print("ploting the composite analysis of surface temperature")
        var_data = self.read_var()

        first_index = self.first_pc.sel(plev=plev)
        last_index = self.last_pc.sel(plev=plev)

        first_var = composite.Tel_field_composite(first_index, var_data,reduction='mean',threshold=1.5)
        last_var = composite.Tel_field_composite(last_index, var_data,reduction='mean',threshold=1.5)

        temp_NAO = composite.composite_plot(first_var, last_var, "NAO")
        plt.savefig(
            self.to_plot_dir + f"composite_tsurf_NAO.png",
            dpi=300,
        )

        temp_EA = composite.composite_plot(first_var, last_var, "EA")
        plt.savefig(
            self.to_plot_dir  + f"composite_tsurf_EA.png",
            dpi=300,
        )

    def plot_all(self):
        self.stat_overview()
        self.spatial_pattern_change()
        self.extreme_count_profile(xlim=(20, 120))
        self.extrc_tsurf(ylim=(25, 130))
        self.extrc_atg(ylim=(25, 130))
        self.composite_analysis()

    def create_doc(self):
        """create md file for the plots"""
        print("creating md files")
        with open(self.doc_dir + self.prefix + "story_line.md", "w") as f:
            f.write(f"# {self.prefix}Story line\n")
            f.write("## Statistical overview\n")
            f.write(
                f"![stat_overview](plots/story_line/{self.prefix}stat_overview.png)\n"
            )
            f.write("## Spatial pattern change at 0K, 2K, 4K\n")
            f.write(
                f"![spatial_pattern_change](plots/story_line/{self.prefix}spatial_pattern_change_500hPa.png)\n"
            )
            f.write("## Extreme event count profile\n")
            f.write(
                f"![extreme_count_vertical_profile](plots/story_line/{self.prefix}extreme_count_vertical_profile.png)\n"
            )
            f.write("## Extreme event count vs. tsurf\n")
            f.write(
                f"![extreme_count_tsurf](plots/story_line/{self.prefix}extreme_count_tsurf_500hPa.png)\n"
            )
            f.write("## Extreme event count vs. arctic-tropical gradient\n")
            f.write(
                f"![extreme_count_atg](plots/story_line/{self.prefix}extreme_count_atg_500hPa.png)\n"
            )
            f.write(
                "## Composite analysis of surface temperature in terms of different extreme events\n"
            )
            f.write("### NAO\n")
            f.write(
                f"![composite_tsurf_NAO](plots/story_line/{self.prefix}composite_tsurf_NAO.png)\n"
            )
            f.write("### EA\n")
            f.write(
                f"![composite_tsurf_EA](plots/story_line/{self.prefix}composite_tsurf_EA.png)\n"
            )
