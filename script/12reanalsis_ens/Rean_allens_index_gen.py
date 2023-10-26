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
def linear_trend(xarr):
    linear_coef = xarr.polyfit(dim="time", deg=1)
    linear_fitted = xr.polyval(xarr.time, linear_coef.polyfit_coefficients)
    return linear_fitted
def quadratic_trend(xarr):
    quadratic_coef = xarr.polyfit(dim="time", deg=2)
    quadratic_fitted = xr.polyval(xarr.time, quadratic_coef.polyfit_coefficients)
    return quadratic_fitted

def detrend(data,method = 'linear_trend'):
    ens_data = data.copy()
    try:
        ens_data = ens_data.mean(dim = 'ens')
    except ValueError:
        pass
    if method == 'linear_trend':
        fitted = ens_data.groupby('time.month').apply(linear_trend)
    elif method == 'quadratic_trend':
        fitted = ens_data.groupby('time.month').apply(quadratic_trend)
    detrended = data - fitted
    return detrended

# %%
# read gph data
def read_gph_data(model, external_forcing="quadratic_trend", **kwargs):
    plev = 50000
    odir = "/work/mh0033/m300883/Tel_MMLE/data/" + model + "/"
    start_year = kwargs.get("start_year", "1940")
    end_year = kwargs.get("end_year", "2022")

    data_JJA = []
    for month in ["Jun", "Jul", "Aug"]:
        print(f"reading the gph data of {month} ...")
        zg_path = odir + "zg_" + month + "/"
        file_names = sorted(glob.glob(zg_path + "*.nc"))
        

        # if model contains 'all_ens'
        if "allens" in model:
            try:
                data_month = xr.open_mfdataset(file_names, combine="nested", concat_dim="ens")
                data_month = data_month['HGT']
            except KeyError:
                data_month = data_month['zg']
        else:
            if model == "CR20":
                data_month = xr.open_dataset(file_names[0])
                data_month = data_month.sel(level = plev/100)
                data_month = data_month['hgt']
            elif model == "ERA5":
                data_month = xr.open_mfdataset(file_names, combine="by_coords")
                data_month = data_month.sel(plev = plev)
                data_month = data_month['var129']
        data_month = data_month.sel(time=slice(start_year, end_year))
        data_JJA.append(data_month)
    data = xr.concat(data_JJA, dim="time").sortby("time")
    # remove the forced trend
    print("detrending ...")
    detrended = detrend(data,method = external_forcing)
    return detrended.squeeze()


# %%
def decompose_period(xarr, nmode = 2):
    xarr = xarr.fillna(0) # fill nan with 0
    field = xarr.sortby("time")
    if 'ens' in xarr.dims:
        field = field.stack(com=("ens", "time"))
        dim = 'com'
    else:
        dim = 'time'
    eof_result = ssp.doeof(field, standard="eof_spatial_std",nmode=nmode,dim = dim)
    return eof_result


# %%
def standard_period(first_eof, last_eof):
    if 'ens' in first_eof.dims:
        mean = first_eof["pc"].mean(dim=("time", "ens"))
        std = first_eof["pc"].std(dim=("time", "ens"))
    else:
        mean = first_eof["pc"].mean(dim="time")
        std = first_eof["pc"].std(dim="time")

    first_eof["pc"] = (first_eof["pc"] - mean) / std
    last_eof["pc"] = (last_eof["pc"] - mean) / std
    return first_eof, last_eof

#%%

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
                f"/work/mh0033/m300883/Tel_MMLE/docs/source/plots/supplyment/{self.model}_{self.group_size}_monthly_spatial_index.png"
            )
#%%
# ERA5, 1950 - 2022, 40 years
ERA5_40 = EOF_reannalyis("ERA5", group_size=40, start_year="1950", end_year="2022")
ERA5_40.save_eof()
ERA5_40.plot_spatial_index(save=True)

#%%
# RA5, 1979 - 2022, 20 years
ERA5_20 = EOF_reannalyis("ERA5", group_size=20, start_year="1979", end_year="2022")
ERA5_20.save_eof()
ERA5_20.plot_spatial_index(save=True)

#%%
# RA5_allens, 1950 - 2022, 40 years
ERA5_allens_40 = EOF_reannalyis("ERA5_allens", group_size=40, start_year="1950", end_year="2022")
ERA5_allens_40.save_eof()
ERA5_allens_40.plot_spatial_index(save=True)
#%%
# CR20, 1950 - 2015, 40 years
CR20_40 = EOF_reannalyis("CR20", group_size=40, start_year="1950", end_year="2015")
CR20_40.save_eof()
CR20_40.plot_spatial_index(save=True)

#%%
# CR20_allens, 1950 - 2015, 40 years
CR20_allens_40 = EOF_reannalyis("CR20_allens", group_size=40, start_year="1950", end_year="2015")
CR20_allens_40.save_eof()
CR20_allens_40.plot_spatial_index(save=True)
# %%
