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
class period_index:
    def __init__(
        self,
        vertical_eof: str,  # 'ind' or 'dep'
        fixed_pattern: str,  # 'all','first','last'
        compare: str,  #'temp', 'CO2'
    ):

        #####################################
        #### some definitions here ##########
        self.vertical_eof = vertical_eof
        self.fixed_pattern = fixed_pattern
        self.compare = compare

        #### some locations here #####
        # the loc for EOF result
        self.eof_dir = (
            "/work/mh0033/m300883/3rdPanel/data/class_decompose/"
            + self.fixed_pattern
            + "Pattern/"
            + self.vertical_eof
            + "/"
        )

        # the loc for the tsurf for determine the time of 1, 2, and 4 degree.
        self.tsurf_dir = "/work/mh0033/m300883/3rdPanel/data/mean_global_temp/mtsurf_mpi_GE_onepct.nc"

        # the destination for savinig plots
        self.plot_dir = (
            "/work/mh0033/m300883/Tele_season/docs/source/plots/quick_plots/"
        )

        # the destination for the doc
        self.img_dir = "plots/quick_plots/"  # relative, no why
        self.doc_dir = "/work/mh0033/m300883/Tele_season/docs/source/"

        ###########################################
        #### tools for naming #####################
        self.prefix = (
            self.vertical_eof + "_" + self.fixed_pattern + "_"
        )  # for name/ ind_all_

        ###############################################
        ##### the data reading and preprocessing ######
        # read data of eof, index and explained variance
        self.eof, self.pc, self.fra = self.read_eof_data()
        self.pc["time"] = self.pc.indexes["time"].to_datetimeindex()

        # data of 500 hpa.
        self.eof_500hpa, self.pc_500hpa, self.fra_500hpa = self.sel_500hpa()
        self.pc_500hpa_df = self.first10_last10_index_df(self.pc_500hpa)

        # read the original gph data to do the composite spatial pattern
        self.gph = self.read_gph_data()

        # index of different period to compare, either first10 v.s last10, or 0,2,4 .C (degree)
        self.pc_periods = self.split_period()

        # index of first10 and last10
        self.first10_pc = self.pc.isel(time=slice(0, 10))
        self.last10_pc = self.pc.isel(time=slice(-10, self.pc.time.size))

        # extreme counts
        self.first_ext_count = extreme.period_extreme_count(self.first10_pc)
        self.last_ext_count = extreme.period_extreme_count(self.last10_pc)

    def read_eof_data(self):
        """
        The data are stored in `3rdPanel/data/`
        """
        print("reading eof result data...")
        odir = self.eof_dir
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
        trop["time"] = trop.indexes[
            "time"
        ].to_datetimeindex()  # the same time type as index
        return trop

    def sel_500hpa(self):
        eof_500 = self.eof.sel(hlayers=50000)
        pc_500 = self.pc.sel(hlayers=50000)

        if self.vertical_eof == "ind":
            fra_500 = self.fra.sel(hlayers=50000)
        elif self.vertical_eof == "dep":
            fra_500 = self.fra

        return eof_500, pc_500, fra_500

    def period_CO2(self):
        """select the first10 and last10 years"""
        first10_pc = self.pc.isel(time=slice(0, 10))
        last10_pc = self.pc.isel(time=slice(-10, self.pc.time.size))
        periods = [first10_pc, last10_pc]
        return periods

    def return_year(self, xarr):
        """return the ten year slice to select"""
        start = xarr.time.values + DateOffset(years=-4)
        end = xarr.time.values + DateOffset(years=5)
        return slice(str(start.year), str(end.year))

    def temp_period(self, fldmean: xr.DataArray):
        """
        to calculate the year when the mean global surface temperature reaches 1,2,and 4 degrees.
        **Argument**
            *fldmean* the fldmean of tsurf
        """
        if isinstance(fldmean, xr.DataArray):
            pass
        else:
            print("only DataArray is accapted, DataSet recevied")

        try:
            fldmean.lon.size == 1 & fldmean.lat.size == 1
        except ValueError:
            print("the fldmean temperature should be calculated first")

        # ens mean
        if fldmean.ens.size != 1:
            fld_ens_mean = fldmean.mean(dim="ens")
        else:
            fld_ens_mean = fldmean

        # squeeze
        mean = fld_ens_mean.squeeze()

        # anomaly
        anomaly = mean - mean[0]

        periods = []

        # 0 degree (1855)
        periods.append(self.return_year(anomaly[5]))

        # 2 degree
        periods.append(
            self.return_year(anomaly.where(anomaly >= 2, drop=True)).squeeze()[0]
        )

        # 4 degree
        periods.append(
            self.return_year(anomaly.where(anomaly >= 4, drop=True)).squeeze()[0]
        )

        return periods

    def CO2_period(self):
        """select the year from pc"""
        first10 = slice("1851", "1860")
        last10 = slice("1990", "1999")
        periods = [first10, last10]
        return periods

    def split_period(self):
        if self.compare == "CO2":
            periods = self.CO2_period()
        elif self.compare == "temp":
            print("reading the mean tsurf data...")
            tsurf = xr.open_dataset(self.tsurf_dir).tsurf
            periods = self.temp_period(tsurf)
        pc_period = []
        for period in periods:
            pc_period.append(self.pc.sel(time=period))
        return pc_period

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


# %%
