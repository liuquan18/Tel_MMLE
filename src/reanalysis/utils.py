import glob
import xarray as xr
import numpy as np


# %%
def linear_trend(xarr):
    linear_coef = xarr.polyfit(dim="time", deg=1)
    linear_fitted = xr.polyval(xarr.time, linear_coef.polyfit_coefficients)
    return linear_fitted


def quadratic_trend(xarr):
    quadratic_coef = xarr.polyfit(dim="time", deg=2)
    quadratic_fitted = xr.polyval(xarr.time, quadratic_coef.polyfit_coefficients)
    return quadratic_fitted


def detrend(data, method="linear_trend"):
    ens_data = data.copy()
    try:
        ens_data = ens_data.mean(dim="ens")
    except ValueError:
        pass
    if method == "linear_trend":
        fitted = ens_data.groupby("time.month").apply(linear_trend)
    elif method == "quadratic_trend":
        fitted = ens_data.groupby("time.month").apply(quadratic_trend)
    elif method == "ens_mean":
        fitted = ens_data.mean(dim="ens")
    detrended = data - fitted
    return detrended


# %%
# read gph data
def read_gph_data(model, external_forcing="quadratic_trend", **kwargs):
    plev = kwargs.get("plev", 50000)
    variable = kwargs.get("variable", "zg")

    odir = "/work/mh0033/m300883/Tel_MMLE/data/" + model + "/"
    start_year = kwargs.get("start_year", "1940")
    end_year = kwargs.get("end_year", "2022")

    data_JJA = []
    for month in ["Jun", "Jul", "Aug"]:
        print(f"reading the gph data of {month} ...")
        zg_path = odir + f"{variable}_" + month + "/"
        file_names = sorted(glob.glob(zg_path + "*.nc"))

        # if model contains 'all_ens'
        if "allens" in model:
            try:
                data_month = xr.open_mfdataset(
                    file_names,
                    combine="nested",
                    concat_dim="ens",
                )
                if variable == 'zg_all_plev':
                    data_month = data_month.load()
                data_month = data_month.sel(plev=plev)
                data_month = data_month["HGT"]
            except KeyError:
                data_month = data_month["zg"]
        else:
            if model == "CR20":
                data_month = xr.open_dataset(file_names[0])
                data_month = data_month.sel(level=plev / 100)
                data_month = data_month["hgt"]
            elif model == "ERA5":
                data_month = xr.open_mfdataset(file_names, combine="by_coords")
                data_month = data_month.sel(plev=plev)
                data_month = data_month["var129"]
        data_month = data_month.sel(time=slice(start_year, end_year))
        data_JJA.append(data_month)
    data = xr.concat(data_JJA, dim="time").sortby("time")
    # remove the forced trend
    if external_forcing is None:
        detrended = data
    else:
        print("detrending ...")
        detrended = detrend(data, method=external_forcing)
    return detrended.squeeze()
