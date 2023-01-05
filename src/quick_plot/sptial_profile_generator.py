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

# functions to plot
import src.plots.vertical_profile as profile_plots
import src.plots.PDF as pdf_plots
import src.plots.plot_violin as violin_plots
import src.plots.spatial_distribution_plot as spatial_dis_plots
import src.plots.return_period as RP_plots
import src.plots.composite_spatial_pattern as composite_spatial_pattern
import src.plots.composite_var as composite_var

import src.extreme.period_pattern_extreme as extreme
import src.EVT.return_period as EVT
import src.composite.field_composite as composite
import src.html.create_md as create_md
import src.Teleconnection.tools as tools


#%%
import importlib

importlib.reload(create_md)

#%%
class first10_last10_index:
    def __init__(
        self,
        vertical_eof: str,  # 'ind' or 'dep'
        fixed_pattern: str,  # 'all','first','last'
    ):

        #### some definitions here ####
        self.vertical_eof = vertical_eof
        self.fixed_pattern = fixed_pattern

        self.prefix = (
            self.vertical_eof + "_" + self.fixed_pattern + "_"
        )  # for name/ ind_all_

        # the destination for savinig plots
        self.plot_dir = "/work/mh0033/m300883/Tele_season/docs/source/plots/quick_plots/"

        # the destination for the doc
        self.img_dir = "plots/quick_plots/"  # relative, no why
        self.doc_dir = "/work/mh0033/m300883/Tele_season/docs/source/"

        # read data of eof, index and explained variance
        self.eof, self.pc, self.fra = self.read_eof_data()
        self.pc["time"] = self.pc.indexes["time"].to_datetimeindex()

        # data of 500 hpa.
        self.eof_500hpa, self.pc_500hpa, self.fra_500hpa = self.sel_500hpa()
        self.pc_500hpa_df = self.first10_last10_index_df(self.pc_500hpa)

        # index of first10 and last10
        self.first10_pc = self.pc.isel(time=slice(0, 10))
        self.last10_pc = self.pc.isel(time=slice(-10, self.pc.time.size))

        # extreme counts
        self.first_ext_count = extreme.period_extreme_count(self.first10_pc)
        self.last_ext_count = extreme.period_extreme_count(self.last10_pc)

        # read the original gph data to do the composite spatial pattern
        self.gph = self.read_gph_data()


    def read_eof_data(self):
        """
        The data are stored in `3rdPanel/data/`
        """
        print("reading eof result data...")
        odir = (
            "/work/mh0033/m300883/3rdPanel/data/class_decompose/"
            + self.fixed_pattern
            + "Pattern/"
            + self.vertical_eof
            + "/"
        )
        eof = xr.open_dataset(odir + self.prefix + "eof.nc").eof

        pc = xr.open_dataset(odir + self.prefix + "pc.nc").pc

        fra = xr.open_dataset(odir + self.prefix + "fra.nc").exp_var
        return eof, pc, fra

    def read_gph_data(self):
        """
        The data are stored somewhere here...
        """
        print("reading gph data...")
        # data
        allens = xr.open_dataset(
            "/work/mh0033/m300883/transition/gr19/gphSeason/allens_season_time.nc"
        )
        # split ens
        splitens = tools.split_ens(allens)

        # demean ens-mean
        demean = splitens - splitens.mean(dim="ens")

        # select traposphere
        trop = demean.sel(hlayers=slice(20000, 100000))

        trop = trop.var156

        trop = tools.standardize(trop)
        trop['time'] = trop.indexes['time'].to_datetimeindex() # the same time type as index
        return trop

    def sel_500hpa(self):
        eof_500 = self.eof.sel(hlayers=50000)
        pc_500 = self.pc.sel(hlayers=50000)

        if self.vertical_eof == "ind":
            fra_500 = self.fra.sel(hlayers=50000)
        elif self.vertical_eof == "dep":
            fra_500 = self.fra

        return eof_500, pc_500, fra_500

    def first10_last10_index_df(self, index):

        """
        select the period data, transform to dataframe
        """
        first10 = index.isel(
            time=slice(0, 10)
        )  # first10 years projectecd on all, standardized with whole time
        last10 = index.isel(time=slice(-10, self.pc.time.size))

        # to dataframe()
        coords = xr.IndexVariable(dims="periods", data=["first10", "last10"])
        index_500hpa = xr.concat([first10, last10], dim=coords)
        index_500hpa = index_500hpa.to_dataframe().reset_index()

        return index_500hpa

    def plot_500hpa_spatial_violin(self):
        """
        sptail maps and violin plots of indexes (NAO and EA).
        """
        print("ploting spatial patterns and violin plot of NAO and EA index ...")
        fig = spatial_dis_plots.spatialMap_violin(
            self.eof_500hpa, self.pc_500hpa_df, self.fra_500hpa
        )

        plt.savefig(
            self.plot_dir + self.prefix + "spatial_pattern_violin500hpa.png", dpi=300
        )

    def plot_500hpa_spatial_hist(self):
        """
        sptail maps and violin plots of indexes.
        """
        print("ploting spatial patterns map and histgram of NAO and EA index ...")
        fig = spatial_dis_plots.spatialMap_hist(
            self.eof_500hpa, self.pc_500hpa_df, self.fra_500hpa
        )

        plt.savefig(
            self.plot_dir + self.prefix + "spatial_pattern_hist500hpa.png", dpi=300
        )

    def violin_profile(self):
        print("ploting the violin profile of NAO and EA index ...")
        fig = violin_plots.plot_vilion(self.first10_pc, self.last10_pc, "whole")
        plt.savefig(self.plot_dir + self.prefix + "violin_profile.png", dpi=300)

    def extreme_count_profile(self, mode):
        print(f"ploting the profile of extreme event count of {mode} index ...")
        fig = profile_plots.plot_vertical_profile(
            self.first_ext_count, self.last_ext_count, mode=mode, std_type="all"
        )
        plt.savefig(
            self.plot_dir + self.prefix + mode + "_extreme_count_profile.png", dpi=300
        )

    def return_period_scatter(self, mode, hlayers=50000):
        print("scatter plot of return period")
        fig = RP_plots.return_period_scatter(self.pc, mode, hlayers=hlayers)
        plt.savefig(
            self.plot_dir + self.prefix + mode + "_return_period_scatter.png", dpi=300
        )

    def return_period_profile(self, mode):
        pos, neg = EVT.vertical_return_period(self.pc, mode)
        fig = RP_plots.return_period_profile(pos, neg, self.pc, mode)
        plt.savefig(
            self.plot_dir + self.prefix + mode + "_return_period_profile.png", dpi=300
        )

    def extreme_spatial_pattern(self, hlayers=100000):
        # do the composite of gph to get the extreme sptial patterns
        first_sptial_pattern = composite.Tel_field_composite(self.first10_pc, self.gph)
        last_sptial_pattern = composite.Tel_field_composite(self.last10_pc, self.gph)
        fig = composite_spatial_pattern.composite_spatial_pattern(
            first_sptial_pattern,
            last_sptial_pattern,
            levels=np.arange(-2, 2.1, 0.4),
            hlayers=hlayers,
        )
        plt.savefig(
            self.plot_dir
            + self.prefix
            + f"extreme_spatial_pattern_{hlayers/100:.0f}hpa.png",
            dpi=300,
        )

    def read_var(self,var):
        """
        read the tsurf or wind, precipitation for composite analysis
        data are stored in 3rdPanel/data/
        """
        data_path = (
        "/work/mh0033/m300883/3rdPanel/data/influence/"
        + var
        + "/"
        + "onepct_1850-1999_ens_1-100."
        + var
        + ".nc"
        )
        var_data = xr.open_dataset(data_path)[var]
        return var_data


    def composite_var(self,var, mode,hlayers = 50000):

        var_data = self.read_var(var)
        var_data = var_data-var_data.mean(dim = 'ens') # demean ens mean

        first_index = self.first10_pc.sel(hlayers = hlayers)
        last_index = self.last10_pc.sel(hlayers = hlayers)

        first_var = composite.Tel_field_composite(first_index,var_data)
        last_var = composite.Tel_field_composite(last_index,var_data)

        fig = composite_var.composite_var(var,first_var, last_var,mode )
        plt.savefig(
            self.plot_dir
            + self.prefix
            + f"composite_{var}_{mode}.png",
            dpi=300,
        )


    def plot_all(self):
        self.plot_500hpa_spatial_violin()
        self.plot_500hpa_spatial_hist()
        self.violin_profile()
        self.extreme_count_profile("NAO")
        self.extreme_count_profile("EA")
        self.return_period_scatter("NAO")
        self.return_period_scatter("EA")
        self.return_period_profile("NAO")
        self.return_period_profile("EA")
        self.extreme_spatial_pattern(hlayers=100000)
        self.composite_var('tsurf','NAO',hlayers=50000)
        self.composite_var('tsurf','EA', hlayers = 50000)

    def create_doc(self):
        create_md.doc_quick_plots(
            self.doc_dir + self.vertical_eof + "_" + self.fixed_pattern,
            f"{self.vertical_eof} decomposition {self.fixed_pattern}-pattern quick plots",
            self.img_dir,
            self.prefix,
        )

if __name__ == "__main__":
    # %%
    ind_all = first10_last10_index("ind", "all")
    ind_all.plot_all()
    ind_all.create_doc()
    # %%
    dep_all = first10_last10_index("dep", "all")
    dep_all.plot_all()
    dep_all.create_doc()
    # %%
