# %%
import xarray as xr
import numpy as np
import os
import sys

import src.MMLE_TEL.index_generator as index_generator
import src.Teleconnection.tools as tools
from scipy.stats import linregress


# %%
def read_gph_data(model,remove_ens_mean = True):
    odir = f"/work/mh0033/m300883/Tel_MMLE/data/{model}/"

    # read gph data
    data_JJA = []
    for month in ["Jun", "Jul", "Aug"]:
        print(f"reading the gph data of {month} ...")
        zg_path = odir + "zg_" + month + "/"
        data_JJA.append(index_generator.read_data(zg_path, plev=50000, remove_ens_mean=remove_ens_mean))
    data = xr.concat(data_JJA, dim="time").sortby("time")
    data = change_lon_to_180(data)
    return data


# function to change the lontitude from 0-360 to -180-180
def change_lon_to_180(zg):
    zg = zg.assign_coords(lon=(((zg.lon + 180) % 360) - 180)).sortby("lon")
    return zg


# %%


def standardize_arr(arr, standard="first10", dim="com"):
    """
    standard when count the extreme cases.
    """

    if standard == "first10":
        years = np.unique(arr.sortby("time").time.dt.year)
        arr_ref = arr.sel(time=slice(str(years[0]), str(years[9])))
        arr_ref = arr_ref.stack(com=("time", "ens"))
        arr_ref_mean = arr_ref.mean(dim=dim)
        arr_ref_std = arr_ref.std(dim=dim)
        arr_standard = (arr - arr_ref_mean) / arr_ref_std
    elif standard == "all":
        arr = arr.stack(com=("time", "ens"))
        arr_box_mean = arr.mean(dim=dim)
        arr_box_std = arr.std(dim=dim)
        arr_standard = (arr - arr_box_mean) / arr_box_std
    return arr_standard


def stats_arr(arr, statis="std", **kwargs):
    if statis == "std":
        arr = arr.stack(com=("time", "ens"))
        arr_stat = arr.std(dim="com")
    elif statis == "mean":
        arr = arr.stack(com=("time", "ens"))
        arr_stat = arr.mean(dim="com")
    elif statis == "count_pos":
        standard = kwargs.get("standard", "first10")
        threshold = kwargs.get("threshold", 1.5)
        arr_standard = standardize_arr(arr, standard=standard)
        arr_standard = arr_standard.stack(com=("time", "ens"))
        # count the number of events above 1.5
        arr_stat = arr_standard.where(arr_standard > threshold).count(dim="com")
    elif statis == "count_neg":
        standard = kwargs.get("standard", "first10")
        threshold = kwargs.get("threshold", 1.5)
        arr_standard = standardize_arr(arr, standard=standard)
        arr_standard = arr_standard.stack(com=("time", "ens"))
        # count the number of events above 1.5
        arr_stat = arr_standard.where(arr_standard < -1 * threshold).count(dim="com")
    elif statis == "slope":
        arr = arr.stack(com=("time", "ens"))
        arr_stat = arr.polyfit(dim="com", deg=1)
    else:
        print("please input the correct statistic method")

    return arr_stat


# %%
def box_spatial_mean(xarr, blat, tlat, llon, rlon):
    """
    calculate the spatial mean of the box.
    check if the lat and lon of the eof is from lower to higher.
    input the order of expected box in # bottom lat, top lat, left lon, right lon for the box
    return the order of the box in the eof.
    """
    if xarr.lat[0] > xarr.lat[-1]:  # descending
        start_lat = tlat
        end_lat = blat
    else:  # ascending
        start_lat = blat
        end_lat = tlat

    if xarr.lon[0] > xarr.lon[-1]:
        start_lon = rlon
        end_lon = llon

    else:
        start_lon = llon
        end_lon = rlon

    arr_box = xarr.sel(lat=slice(start_lat, end_lat), lon=slice(start_lon, end_lon))
    weights = np.cos(np.deg2rad(arr_box.lat))
    region_box_weighted = arr_box.weighted(weights)
    arr_box_spatial_mean = region_box_weighted.mean(("lon", "lat"))

    return arr_box_spatial_mean


def NAO_pos_center(zg, statis="std"):
    pos_box_spa_mean = box_spatial_mean(zg, 45, 60, -30, 0)
    pos_box_std = stats_arr(pos_box_spa_mean, statis=statis)
    return pos_box_std


def NAO_neg_center(zg, statis="std"):
    neg_box = box_spatial_mean(zg, 60, 75, -75, -50)
    neg_box_std = stats_arr(neg_box, statis=statis)
    return neg_box_std
    # spaital mean with weights


def box_variability(model):
    zg = read_gph_data(model)

    pos_var = zg.resample(time="10AS-JUN").apply(NAO_pos_center, statis="std")
    neg_var = zg.resample(time="10AS-JUN").apply(NAO_neg_center, statis="std")
    return pos_var, neg_var


def box_mean(model):
    zg = read_gph_data(model)

    pos_var = zg.resample(time="10AS-JUN").apply(NAO_pos_center, statis="mean")
    neg_var = zg.resample(time="10AS-JUN").apply(NAO_neg_center, statis="mean")
    return pos_var, neg_var


# %%
# define function to calculate slope and p-value
def linregress_ufunc(x):
    slope, _, r_value, p_value, _ = linregress(np.arange(len(x)), x)
    return slope, p_value < 0.05


# %%
def slope_ens_std(model):
    zg = read_gph_data(model)
    zg.load()
    zg_std = zg.resample(time="10AS-JUN").apply(stats_arr, statis="std")

    result = xr.apply_ufunc(
        linregress_ufunc,
        zg_std,
        input_core_dims=[["time"]],
        output_core_dims=[[], []],
        vectorize=True,
    )

    # convert result to DataArray
    slope_da = result[0].rename("slope")
    pvalue_da = result[1].rename("pvalue")

    # create a new dataset and add the slope_da and pvalue_da as variables
    result_ds = xr.Dataset({"slope": slope_da, "pvalue": pvalue_da})
    return result_ds
# %%

def slope_ens_mean(model):
    zg = read_gph_data(model,remove_ens_mean=False)
    zg.load()

    # calculate the decadal mean of zg to corresponds to the slope of ens std.
    zg_yearly = zg.resample(time="10AS-JUN").mean(dim="time")

    # calculate the slope and pvalue
    result = xr.apply_ufunc(
        linregress_ufunc,
        zg_yearly,
        input_core_dims=[["time"]],
        output_core_dims=[[], []],
        vectorize=True,
    )

    # convert result to DataArray
    slope_da = result[0].rename("slope")
    pvalue_da = result[1].rename("pvalue")

    # create a new dataset and add the slope_da and pvalue_da as variables
    result_ds = xr.Dataset({"slope": slope_da, "pvalue": pvalue_da})
    return result_ds

#%%