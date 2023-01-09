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
        model:str, # the model name
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

        #### some locations here #####
        odir = "/work/mh0033/m300883/Tel_MMLE/data/"+self.model
        # the loc for EOF result
        self.eof_dir = odir+"/EOF_result/"

        # the loc for the tsurf for determine the time of 1, 2, and 4 degree.
        self.tsurf_fldmean_dir = odir + "/ts_preprocessed/"

        # the loc for original geopotential height data
        self.zg_dir = odir + "/zg/"

        # the loc for original tsurface map
        self.ts_dir = odir + "/ts/"

        # the destination for savinig plots
        self.plot_dir = (
            "/work/mh0033/m300883/Tel_MMLE/docs/source/plots/" + self.model + "/"
        )

        # the destination for the doc
        self.img_dir = "plots/" + self.model + "/"  # relative, no why
        self.doc_dir = "/work/mh0033/m300883/Tel_MMLE/docs/source/"

        ###########################################
        #### tools for naming #####################
        self.prefix = (
            self.vertical_eof + "_" + self.fixed_pattern + "_"
        )  # for name/ ind_all_
        if self.compare == "CO2":
            self.period_name = ["first10", "last10"]
        elif self.compare == "temp":
            self.period_name = ["0C", "2C", "4C"]

        ###############################################
        ##### the data reading and preprocessing ######
        # read data of eof, index and explained variance
        self.eof, self.pc, self.fra = self.read_eof_data()
        self.pc["time"] = self.pc.indexes["time"].to_datetimeindex()

        # read the original gph data to do the composite spatial pattern
        self.gph = self.read_gph_data()

        # index of different period to compare, either first10 v.s last10, or 0,2,4 .C (degree)
        self.pc_periods, self.periods = self.split_period()

        # extreme counts
        self.ext_counts_periods = self.extreme_periods()

        # data of 500 hpa.
        self.eof_500hpa, self.pc_500hpa, self.fra_500hpa = self.sel_500hpa()
        self.pc_500hpa_df = self.bar500hpa_index_df()

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
        zg_data = xr.open_mfdataset(self.zg_dir,combine="nested",concat_dim="ens")

        # demean ens-mean
        demean = zg_data - zg_data.mean(dim="ens")

        # select traposphere
        if self.mode == "MPI_GE":
            trop = demean.sel(hlayers=slice(20000, 100000))
            trop = trop.var156
        else:
            trop = demean.sel(hlayers=slice(100000,20000))
            trop = trop.zg

        trop = tools.standardize(trop)
        trop["time"] = trop.indexes["time"].to_datetimeindex() 
        return trop

    def read_var(self, var):
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
            self.return_year(anomaly.where(anomaly >= 2, drop=True).squeeze()[0])
        )

        # 4 degree
        periods.append(
            self.return_year(anomaly.where(anomaly >= 4, drop=True).squeeze()[0])
        )

        return periods

    def CO2_period(self):
        """select the year from pc"""
        years = self.pc.time
        first10 = slice(years[0], years[10])
        last10 = slice(years[-10], years[years.size])
        periods = [first10, last10]
        return periods

    def split_period(self):
        if self.compare == "CO2":
            periods = self.CO2_period()
        elif self.compare == "temp":
            print("reading the mean tsurf data...")
            if self.model == "MPI_GE":
                tsurf = xr.open_dataset(self.tsurf_fldmean_dir)  # already pre-processed
            else:
                tsurf = xr.open_mfdataset(self.tsurf_fldmean_dir,combine = 'nested',concat_dim = 'ens')
            try: # different name
                tsurf = tsurf.tsurf
            except AttributeError:
                tsurf = tsurf.ts
            periods = self.temp_period(tsurf)
        pcs_period = []
        for i, period in enumerate(periods):
            pc_period = self.pc.sel(time=period)
            pc_period["compare"] = self.period_name[i]
            pcs_period.append(pc_period)
        return pcs_period, periods

    def extreme_periods(self):
        ext_counts_list = []
        for pc_period in self.pc_periods:
            ext_counts_list.append(extreme.period_extreme_count(pc_period))
            ext_counts = xr.concat(ext_counts_list, dim="compare")
        return ext_counts

    def bar500hpa_index_df(self):

        """
        select the period data, transform to dataframe
        """
        first = self.pc_periods[0].sel(hlayers=50000)
        last = self.pc_periods[-1].sel(hlayers=50000)

        # to dataframe()
        if self.compare == "CO2":
            coords = xr.IndexVariable(dims="periods", data=["first10", "last10"])
        else:
            coords = xr.IndexVariable(dims="periods", data=["0C", "4C"])
        index_500hpa = xr.concat([first, last], dim=coords)
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
        fig = violin_plots.plot_vilion(
            self.pc_periods[0], self.pc_periods[-1], compare=self.compare
        )
        plt.savefig(self.plot_dir + self.prefix + "violin_profile.png", dpi=300)

    def extreme_count_profile(self, mode):
        print(f"ploting the profile of extreme event count of {mode} index ...")
        fig = profile_plots.plot_vertical_profile(
            self.ext_counts_periods,
            mode=mode,
        )
        plt.savefig(
            self.plot_dir + self.prefix + mode + "_extreme_count_profile.png", dpi=300
        )

    def return_period_scatter(self, mode, hlayers=50000):
        print("scatter plot of return period")
        pos, mpos, neg, mneg = EVT.mode_return_period(
            self.pc, mode=mode, periods=self.periods, hlayers=hlayers
        )
        fig = RP_plots.return_period_scatter(
            pos, mpos, neg, mneg, mode, self.periods,self.period_name, hlayers=hlayers
        )
        plt.savefig(
            self.plot_dir + self.prefix + mode + "_return_period_scatter.png", dpi=300
        )

    def return_period_profile(self, mode):
        pos, neg = EVT.vertical_return_period(self.pc, mode, self.periods)
        fig = RP_plots.return_period_profile(pos, neg, self.pc, mode,self.period_name)
        plt.savefig(
            self.plot_dir + self.prefix + mode + "_return_period_profile.png", dpi=300
        )

    def extreme_spatial_pattern(self, hlayers=100000):
        # do the composite of gph to get the extreme sptial patterns
        first_sptial_pattern = composite.Tel_field_composite(self.pc_periods[0], self.gph)
        last_sptial_pattern = composite.Tel_field_composite(self.pc_periods[-1], self.gph)
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


    def composite_var(self, var, mode, hlayers=50000):

        var_data = self.read_var(var)
        var_data = var_data - var_data.mean(dim="ens")  # demean ens mean

        first_index = self.pc_periods[0].sel(hlayers=hlayers)
        last_index = self.pc_periods[-1].sel(hlayers=hlayers)

        first_var = composite.Tel_field_composite(first_index, var_data)
        last_var = composite.Tel_field_composite(last_index, var_data)

        fig = composite_var.composite_var(var, first_var, last_var, mode)
        plt.savefig(
            self.plot_dir + self.prefix + f"composite_{var}_{mode}.png",
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
        self.composite_var("tsurf", "NAO", hlayers=50000)
        self.composite_var("tsurf", "EA", hlayers=50000)

    def create_doc(self):
        create_md.doc_quick_plots(
            self.doc_dir + self.vertical_eof + "_" + self.fixed_pattern,
            f"{self.vertical_eof} decomposition {self.fixed_pattern}-pattern quick plots",
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