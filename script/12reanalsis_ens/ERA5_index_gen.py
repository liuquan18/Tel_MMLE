# %%
import xarray as xr
import numpy as np
import pandas as pd
import random
import os

import src.warming_stage.warming_stage as warming_stage
import src.Teleconnection.spatial_pattern as ssp
import src.MMLE_TEL.index_generator as index_generate

import glob

import proplot as pplt
import src.plots.extreme_plot as extplt
import src.plots.statistical_overview as stat_overview
import matplotlib.pyplot as plt
import src.Teleconnection.spatial_pattern as ssp
# %%
import importlib

importlib.reload(index_generate)
importlib.reload(ssp)

# %%
def linear_detrend(data):
    ens_mean = data.mean(dim="ens")
    ens_mean_c = ens_mean.copy()
    ens_mean_c["time"] = np.arange(ens_mean_c.time.size)
    linear_coef = ens_mean_c.polyfit(dim="time", deg=1)
    linear_fitted = xr.polyval(ens_mean_c.time, linear_coef.polyfit_coefficients)
    linear_fitted["time"] = ens_mean.time
    linear_detrend = data - linear_fitted
    return linear_detrend

def sel_20CR(xarray):
    return xarray.sel(time=slice("1836", "2015"))
# %%
# read gph data
def read_gph_data(model, external_forcing="linear_trend", **kwargs):
    plev = 50000
    odir = "/work/mh0033/m300883/Tel_MMLE/data/" + model + "/"
    start_year = kwargs.get("start_year", "1940")
    end_year = kwargs.get("end_year", "2022")

    data_JJA = []
    for month in ["Jun", "Jul", "Aug"]:
        print(f"reading the gph data of {month} ...")
        zg_path = odir + "zg_" + month + "/"
        try:
            data_month = index_generate.read_data(zg_path, plev=plev, remove_ens_mean=False)
        except ValueError:
            data_month = xr.open_mfdataset(zg_path + "*.nc", combine="nested", concat_dim="ens",preprocess=sel_20CR)
            data_month = data_month['HGT']
        data_month = data_month.sel(time=slice(start_year, end_year))
        if external_forcing == "linear_trend":
            data_internal = linear_detrend(data_month)
        elif external_forcing == "ens_mean":
            data_internal = data_month - data_month.mean(dim="ens")
        data_JJA.append(data_internal)
    data = xr.concat(data_JJA, dim="time").sortby("time")
    return data


# %%
def decompose_period(xarr, nmode = 2):
    xarr = xarr.fillna(0) # fill nan with 0
    field = xarr.sortby("time")
    field = field.stack(com=("ens", "time"))
    eof_result = ssp.doeof(field, standard="eof_spatial_std",nmode=nmode)
    return eof_result


# %%
def standard_period(first_40_eof, last_40_eof):
    mean = first_40_eof["pc"].mean(dim=("time", "ens"))
    std = first_40_eof["pc"].std(dim=("time", "ens"))

    first_40_eof["pc"] = (first_40_eof["pc"] - mean) / std
    last_40_eof["pc"] = (last_40_eof["pc"] - mean) / std
    return first_40_eof, last_40_eof


class EOF_reannalyis:
    def __init__(self, model, group_size=40, external_forcing="linear_trend", **kwargs):
        self.model = model
        self.group_size = group_size
        self.external_forcing = external_forcing
        self.plev = 50000
        self.standard = "first"
        self.nmode = kwargs.get("nmode", 2)

        self.start_year = kwargs.get("start_year", "1940")
        self.end_year = kwargs.get("end_year", "2022")

        self.data = read_gph_data(model, external_forcing=external_forcing,
                                  start_year = self.start_year, end_year = self.end_year)


        self.first_data = self.data.sel(
            time=slice(self.start_year, str(int(self.start_year) + self.group_size - 1))
        )
        self.last_data = self.data.sel(
            time=slice(str(int(self.end_year) - self.group_size + 1), self.end_year)
        )

        print("********* decomposing *********")
        self.first_eof = decompose_period(self.first_data, nmode=self.nmode)
        self.last_eof = decompose_period(self.last_data, nmode=self.nmode)

        print("********* standardizing *********")
        self.first_eof_std, self.last_eof_std = standard_period(
            self.first_eof, self.last_eof
        )

    def save_eof(self):
        print("********* saving *********")
        self.first_eof.to_netcdf(
            f"/work/mh0033/m300883/Tel_MMLE/data/{self.model}/EOF_result/first_{self.group_size}_eof.nc"
        )
        self.last_eof.to_netcdf(
            f"/work/mh0033/m300883/Tel_MMLE/data/{self.model}/EOF_result/last_{self.group_size}_eof.nc"
        )

        self.first_eof_std.to_netcdf(
            f"/work/mh0033/m300883/Tel_MMLE/data/{self.model}/EOF_result/first_{self.group_size}_eof_std.nc"
        )
        self.last_eof_std.to_netcdf(
            f"/work/mh0033/m300883/Tel_MMLE/data/{self.model}/EOF_result/last_{self.group_size}_eof_std.nc"
        )

    def plot_spatial_index(self,save=False):
        
        fig2 = pplt.figure(figsize=(180 / 25.4, 90 / 25.4), sharex=False, sharey=False)
        fig2.format(
            abc=True,
            abcloc="ul",
            abcstyle="a",
        )
        gs = pplt.GridSpec(
            ncols=2,
            nrows=1,
        )
        ax1 = fig2.add_subplot(gs[0], proj="ortho", proj_kw={"lon_0": -20, "lat_0": 60})
        ax2 = fig2.add_subplot(gs[1])


        ax1, fmap, lmap = stat_overview.spatial_pattern_plot(
            ax1,
            self.first_eof_std.eof.sel(mode="NAO").squeeze(),
            self.first_eof_std.fra.sel(mode="NAO").squeeze(),
            self.last_eof_std.eof.sel(mode="NAO", ).squeeze(),
            self.last_eof_std.fra.sel(mode="NAO").squeeze(),
            levels=np.arange(-2, 2.1, 0.4),
        )
        ax2, hist = stat_overview.index_distribution_plot(
            ax2,
            self.first_eof_std.pc.sel(mode="NAO"),
            self.last_eof_std.pc.sel(mode="NAO"),
        )

        if save:
            plt.savefig(
                f"/work/mh0033/m300883/Tel_MMLE/docs/source/plots/supplyment/{self.model}_monthly_spatial_index.png"
            )

# %%
ERA5_40 = EOF_reannalyis("ERA5_allens", group_size=40, start_year="1940", end_year="2022")
ERA5_40.plot_spatial_index(save=False)
# %%
ERA5_20 = EOF_reannalyis("ERA5_allens", group_size=20, start_year="1940", end_year="2022",nmode = 2)
ERA5_20.plot_spatial_index(save=False)
# %%
ERA5_30 = EOF_reannalyis("ERA5_allens", group_size=30, start_year="1940", end_year="2022",nmode = 2)
ERA5_30.plot_spatial_index(save=False)
# %%
CR20_40 = EOF_reannalyis("CR20", group_size=40, start_year="1930", end_year="2012",nmode = 2)
# %%
CR20_40.plot_spatial_index(save=False)
# %%
CR20_40_allt = EOF_reannalyis("CR20", group_size=40, start_year="1850", end_year="2010",nmode = 2)
CR20_40_allt.plot_spatial_index(save=True)

# %%
CR20_10_allt = EOF_reannalyis("CR20", group_size=10, start_year="1850", end_year="2010",nmode = 2)
CR20_10_allt.plot_spatial_index(save=True)
# %%
