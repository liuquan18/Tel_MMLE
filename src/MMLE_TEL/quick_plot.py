"""
This script generat the plots for:
- the spatial patterns and distribution of index at 500 hpa
- the violin plots
- the vertical profile of extreme counts
- the return period of 500hpa
- the vertical profile of media return period
"""
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
class period_index:
    def __init__(
        self,
        model: str,  # the model name
        vertical_eof: str,  # 'ind' or 'dep'
        fixed_pattern: str,  # 'all','first','last'
        compare: str,  #'temp', 'CO2'
    ):

        #####################################
        #### some definitions here ##########
        self.vertical_eof = vertical_eof
        self.fixed_pattern = fixed_pattern
        self.compare = compare
        self.model = model
        self.prefix = self.vertical_eof + "_" + self.fixed_pattern + "_"

        #### some locations here #####
        odir = "/work/mh0033/m300883/Tel_MMLE/data/" + self.model + "/"
        # the loc for EOF result
        self.eof_dir = odir + "EOF_result/"

        # the loc for the tsurf for determine the time of 1, 2, and 4 degree.
        self.tsurf_fldmean_dir = odir + "ts_processed/"

        # the loc for original geopotential height data
        self.zg_dir = odir + "zg/"

        # the loc for original tsurface map
        self.ts_dir = odir + "ts/"

        # the loc for zg_processed
        self.zg_processed_dir = odir + "zg_processed/"

        # the loc of eofs and fras at all periods
        self.all_eofs_dir = odir + "EOF_result/" + self.prefix + "eofs_allPeriods.nc"
        self.all_fras_dir = odir + "EOF_result/" + self.prefix + "fras_allPeriods.nc"

        # the destination for savinig plots
        self.plot_dir = (
            "/work/mh0033/m300883/Tel_MMLE/docs/source/plots/"
            + self.model
            + "/"
            + self.model
            + "_"
        )

        # the destination for the doc
        self.img_dir = (
            "plots/" + self.model + "/" + self.model + "_"
        )  # relative, no why
        self.doc_dir = "/work/mh0033/m300883/Tel_MMLE/docs/source/"

        ###########################################
        #### tools for naming #####################
        if self.compare == "CO2":
            self.period_name = ["first10", "last10"]
        elif self.compare == "temp":
            self.period_name = ["0K", "2K", "4K"]

        ###############################################
        ##### the data reading and preprocessing ######
        # read data of eof, index and explained variance
        eof_result = self.read_eof_data()
        self.eof = eof_result.eof
        self.fra = eof_result.fra
        self.pc = eof_result.pc

        # fldmean tsurf
        self.fldmean_tsurf = warming_stage.read_tsurf_fldmean(
            self.tsurf_fldmean_dir + "tsurf_mean.nc"
        )

        # index of different period to compare, either first10 v.s last10, or 0,2,4 .K (degree)
        self.pc_periods, self.periods = warming_stage.split_period(
            self.pc, self.compare, self.fldmean_tsurf
        )
        # extreme counts
        self.ext_counts_periods = self.period_wise_extreme()

    def period_wise_extreme(self):
        ext_counts_list = []
        for pc_period in self.pc_periods:
            ext_counts_list.append(extreme.period_extreme_count(pc_period))
            ext_counts = xr.concat(ext_counts_list, dim="compare")
        return ext_counts

    def read_eof_data(self):
        """
        The data are stored in `3rdPanel/data/`
        """
        print("reading eof result data...")
        odir = self.eof_dir
        eof_result = xr.open_dataset(odir + self.prefix + "eof_result.nc")
        return eof_result

    def read_gph_data(self):
        """
        The data are stored somewhere here...
        """
        print("reading gph data...")
        # data
        if self.model == "MPI_GE_onepct":
            zg_data = xr.open_dataset(self.zg_dir + "allens_zg.nc")
            zg_data = tools.split_ens(zg_data)  # not splited yet here

        else:
            zg_data = xr.open_mfdataset(
                self.zg_dir + "*.nc", combine="nested", concat_dim="ens"
            )
            zg_data = zg_data.rename({"plev": "plev"})  # historical error

        # demean ens-mean
        demean = zg_data - zg_data.mean(dim="ens")

        # select traposphere
        trop = demean.sel(plev=slice(100000, 20000))

        if self.model in "MPI_GE_onepct":
            trop = trop.var156
        else:
            trop = trop.zg

        trop = tools.standardize(trop)
        try:
            trop["time"] = trop.indexes["time"].to_datetimeindex()
        except AttributeError:
            trop["time"] = pd.to_datetime(trop.time)
        return trop

    def read_var(self, var):
        """
        read the tsurf or wind, precipitation for composite analysis
        data are stored in 3rdPanel/data/
        """
        print("reading the var to do composite...")
        # MPI is different format...
        if self.model == "MPI_GE_onepct":
            var_data = xr.open_dataset(self.ts_dir + "all_ens_tsurf.nc").tsurf
        elif self.model == "MPI_GE":
            var_data = xr.open_mfdataset(
                self.ts_dir + "*.nc", combine="nested", concat_dim="ens"
            ).tsurf
        else:
            var_data = xr.open_mfdataset(
                self.ts_dir + "*.nc", combine="nested", concat_dim="ens"
            )
            var_data = var_data.ts
            var_data = var_data.rename({"plev": "plev"})
            var_data["time"] = var_data.indexes["time"].to_datetimeindex()

        return var_data

    def sel_500hpa(self):
        eof_500 = self.eof.sel(plev=50000)
        pc_500 = self.pc.sel(plev=50000)

        if self.vertical_eof == "ind":
            fra_500 = self.fra.sel(plev=50000)
        elif self.vertical_eof == "dep":
            fra_500 = self.fra

        return eof_500, pc_500, fra_500

    def bar500hpa_index_df(self):

        """
        select the period data, transform to dataframe
        """
        first = self.pc_periods[0].sel(plev=50000)
        last = self.pc_periods[-1].sel(plev=50000)

        # to dataframe()
        if self.compare == "CO2":
            coords = xr.IndexVariable(dims="periods", data=["first10", "last10"])
        else:
            coords = xr.IndexVariable(dims="periods", data=["0K", "4K"])
        index_500hpa = xr.concat([first, last], dim=coords)
        index_500hpa = index_500hpa.to_dataframe().reset_index()

        return index_500hpa

    #%%
    def plot_500hpa_spatial_violin(self):
        """
        sptail maps and violin plots of indexes (NAO and EA).
        """
        print("ploting spatial patterns and violin plot of NAO and EA index ...")

        # data of 500 hpa.
        eof_500hpa, _, fra_500hpa = self.sel_500hpa()
        try:
            eof_500hpa = eof_500hpa.isel(decade=1)
            fra_500hpa = fra_500hpa.isel(decade=1)
        except ValueError:
            pass

        pc_500hpa_df = self.bar500hpa_index_df()

        fig = spatial_dis_plots.spatialMap_violin(eof_500hpa, pc_500hpa_df, fra_500hpa)

        plt.savefig(
            self.plot_dir + self.prefix + "spatial_pattern_violin500hpa.png", dpi=300
        )

    def plot_500hpa_spatial_hist(self):
        """
        sptail maps and violin plots of indexes.
        """
        print("ploting spatial patterns map and histgram of NAO and EA index ...")

        # data of 500 hpa.
        eof_500hpa, _, fra_500hpa = self.sel_500hpa()
        try:
            eof_500hpa = eof_500hpa.isel(decade=1)
            fra_500hpa = fra_500hpa.isel(decade=1)
        except ValueError:
            pass
        pc_500hpa_df = self.bar500hpa_index_df()

        fig = spatial_dis_plots.spatialMap_hist(eof_500hpa, pc_500hpa_df, fra_500hpa)

        plt.savefig(
            self.plot_dir + self.prefix + "spatial_pattern_hist500hpa.png", dpi=300
        )

    def violin_profile(self):
        print("ploting the violin profile of NAO and EA index ...")
        fig = violin_plots.plot_vilion(
            self.pc_periods[0], self.pc_periods[-1], compare=self.compare
        )
        plt.savefig(self.plot_dir + self.prefix + "violin_profile.png", dpi=300)

    def spatial_change(self):
        print("decomposing spatial patterns at all periods...")
        data = sp_change.read_gph_data(self.zg_processed_dir)
        EOFs, FRAs = sp_change.spatial_pattern_change(
            data, periods=self.periods, names=self.period_name
        )
        return EOFs, FRAs

    def spatial_pattern_change(self):
        print("ploting the spatial pattern changes...")

        # try read data first
        try:
            EOFs = xr.open_dataset(self.all_eofs_dir).eof
            FRAs = xr.open_dataset(self.all_fras_dir).exp_var
            print("found the EOFs and FRAs files")

        except FileNotFoundError:
            print("EOFs file did found, starting a long decomposition now")
            EOFs, FRAs = self.spatial_change()

        maps = sp_change.spatial_pattern_maps(
            EOFs, FRAs, levels=np.arange(-1, 1.1, 0.2)
        )
        plt.savefig(
            self.plot_dir + self.prefix + "spatial_pattern_change_map.png", dpi=300
        )

        # vetmaps = sp_change.spatial_pattern_profile(
        #     EOFs, levels=np.arange(-1.0, 1.1, 0.2)
        # )
        # plt.savefig(
        #     self.plot_dir + self.prefix + "spatial_pattern_change_profile.png", dpi=300
        # )

    def extreme_count_profile(self, mode):
        print(f"ploting the profile of extreme event count of {mode} index ...")
        if self.model not in "MPI_GE_onepct":
            xlim = (0, 20)
        else:
            xlim = (-5, 45)
        fig = profile_plots.plot_vertical_profile(
            self.ext_counts_periods, mode=mode, xlim=xlim
        )
        plt.savefig(
            self.plot_dir + self.prefix + mode + "_extreme_count_profile.png", dpi=300
        )

    def extrc_tsurf_scatter(self, average=False):
        """
        scatter plot of extreme_count v.s fldmean tsurf
        """
        extr_count, ts_mean = extrc_tsurf.decadal_extrc_tsurf(
            self.pc, self.fldmean_tsurf, plev=50000
        )
        if average:
            extr_count = extr_count / self.pc.ens.size

        fig = extrc_tsurf.extCount_tsurf_scatter(extr_count, ts_mean)
        plt.savefig(
            self.plot_dir + self.prefix + "extrc_fldmean_ts_scatter.png", dpi=300
        )

    def return_period_scatter(self, mode, plev=50000):
        print("scatter plot of return period")
        pos, mpos, neg, mneg = EVT.mode_return_period(
            self.pc, mode=mode, periods=self.periods, plev=plev
        )
        fig = RP_plots.return_period_scatter(
            pos, mpos, neg, mneg, mode, self.periods, self.period_name, plev=plev
        )
        plt.savefig(
            self.plot_dir + self.prefix + mode + "_return_period_scatter.png", dpi=300
        )

    def return_period_profile(self, mode):
        pos, neg = EVT.vertical_return_period(self.pc, mode, self.periods)
        fig = RP_plots.return_period_profile(pos, neg, self.pc, mode, self.period_name)
        plt.savefig(
            self.plot_dir + self.prefix + mode + "_return_period_profile.png", dpi=300
        )

    #%%
    def extreme_spatial_pattern(self, plev=100000):

        # read the original gph data to do the composite spatial pattern
        gph = self.read_gph_data()

        # do the composite of gph to get the extreme sptial patterns
        first_sptial_pattern = composite.Tel_field_composite(self.pc_periods[0], gph)
        last_sptial_pattern = composite.Tel_field_composite(self.pc_periods[-1], gph)
        fig = composite_spatial_pattern.composite_spatial_pattern(
            first_sptial_pattern,
            last_sptial_pattern,
            levels=np.arange(-2, 2.1, 0.4),
            plev=plev,
        )
        plt.savefig(
            self.plot_dir
            + self.prefix
            + f"extreme_spatial_pattern_{plev/100:.0f}hpa.png",
            dpi=300,
        )

    def composite_var(self, var, mode, plev=50000):

        var_data = self.read_var(var)
        var_data = var_data - var_data.mean(dim="ens")  # demean ens mean

        first_index = self.pc_periods[0].sel(plev=plev)
        last_index = self.pc_periods[-1].sel(plev=plev)

        first_var = composite.Tel_field_composite(first_index, var_data)
        last_var = composite.Tel_field_composite(last_index, var_data)

        fig = composite_plot.composite_var(var, first_var, last_var, mode)
        plt.savefig(
            self.plot_dir + self.prefix + f"composite_{var}_{mode}.png",
            dpi=300,
        )

    def plot_all(self):
        self.plot_500hpa_spatial_violin()
        self.plot_500hpa_spatial_hist()
        self.violin_profile()
        self.spatial_pattern_change()
        self.extreme_count_profile("NAO")
        self.extreme_count_profile("EA")
        self.extrc_tsurf_scatter()
        self.return_period_scatter("NAO")
        self.return_period_scatter("EA")
        self.return_period_profile("NAO")
        self.return_period_profile("EA")
        # self.extreme_spatial_pattern(plev=100000)
        # self.composite_var("tsurf", "NAO", plev=50000)
        # self.composite_var("tsurf", "EA", plev=50000)

    def create_doc(self):
        create_md.doc_quick_plots(
            self.doc_dir + self.model + "_" + self.prefix + "doc",
            f"{self.model} {self.vertical_eof} decomposition {self.fixed_pattern}-pattern quick plots",
            self.img_dir,
            self.prefix,
        )


#%%
if __name__ == "__main__":
    ind_all = period_index("ind", "all", "CO2")
    ind_all.plot_all()
    ind_all.create_doc()
    dep_all = period_index("dep", "all", "CO2")
    dep_all.plot_all()
    dep_all.create_doc()

# %%
