#%%
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

import src.plots.month_extrc_dig as month_extrc_dig
from scipy.stats import linregress


#%%
import importlib
importlib.reload(month_extrc_dig)

# %%
def read_extrc_month(month, model="MPI_GE_onepct"):
    dir = f"/work/mh0033/m300883/Tel_MMLE/data/{model}/extreme_count/"
    filename = f"plev_50000_decade_first_{month}_extre_counts.nc"
    ds = xr.open_dataset(dir + filename).sel(confidence="true")
    ds["time"] = ds["time.year"]
    return ds.pc


def read_extrc(model="MPI_GE_onepct"):
    months = [
        "Dec", # start from Dec
        "Jan",
        "Feb",
        "Mar",
        "Apr",
        "May",
        "Jun",
        "Jul",
        "Aug",
        "Sep",
        "Oct",
        "Nov",
    ]
    mon_idx = xr.IndexVariable("month", months)
    ds = [read_extrc_month(month, model) for month in months]
    ds = xr.concat(ds, dim=mon_idx)
    return ds

def read_tsurf(model="MPI_GE_onepct"):
    dir = f"/work/mh0033/m300883/Tel_MMLE/data/{model}/extreme_count/"
    filename = f"ens_fld_year_mean.nc"
    ds = xr.open_dataset(dir + filename)
    # ds drop all the dimensions except time
    ds = ds.reset_coords(drop=[coord for coord in ds.coords if coord != "time"])
    ds["time"] = ds["time.year"]
    try:
        ds = ds.tsurf
    except AttributeError:
        ds = ds.ts
    return ds

# %%
# MMLEA
models = ["MPI_GE", "CanESM2", "CESM1_CAM5", "MK36", "GFDL_CM3"]
models_ind = xr.IndexVariable("model", models)

#%%
MMLEA = [read_extrc(model) for model in models]
MMLEA = xr.concat(MMLEA, dim=models_ind)
MMLEA = MMLEA.sel(time=slice("1950-01-01", "2091-01-01"))
#%%
# Calculate the slope of the linear regression along the 'time' dimension
MMLEA_slope = MMLEA.polyfit(dim="time", deg=1).polyfit_coefficients.sel(degree=1)

#%%
Tsurf = [read_tsurf(model) for model in models]
Tsurf = xr.concat(Tsurf, dim=models_ind)
Tsurf = Tsurf.sel(time=slice("1950-01-01", "2091-01-01"))

#%%
def reg_slope(tsurf, extrc):
    slope, _, _, _, _ = linregress(tsurf, extrc)
    return slope

Regress_slope = xr.apply_ufunc(
    reg_slope,
    Tsurf,
    MMLEA,
    input_core_dims=[["time"], ["time"]],
    output_core_dims=[[]],
    vectorize=True,
    dask="parallelized",
    output_dtypes=[float],
)


# %%
importlib.reload(month_extrc_dig)

month_extrc_dig.bar_hatch(MMLEA_slope)

#%%
plt.savefig(
    "/work/mh0033/m300883/Tel_MMLE/docs/source/plots/monthly/month_extrc_diagram_patch_time.png"
)
# %%
month_extrc_dig.bar_side(MMLEA_slope)
plt.savefig(
    "/work/mh0033/m300883/Tel_MMLE/docs/source/plots/monthly/month_extrc_diagram_patch_side_time.png"
)
# %%
month_extrc_dig.bar_hatch(Regress_slope,xlim = (-10,18))
plt.savefig(
    "/work/mh0033/m300883/Tel_MMLE/docs/source/plots/monthly/month_extrc_diagram_patch_tsurf.png"
)
# %%
month_extrc_dig.bar_side(Regress_slope,xlim = (-5,15))
plt.savefig(
    "/work/mh0033/m300883/Tel_MMLE/docs/source/plots/monthly/month_extrc_diagram_patch_side_tsurf.png"
)
# %%
importlib.reload(month_extrc_dig)

#%%
month_extrc_dig.bar_dual(MMLEA_slope)
plt.savefig(
    "/work/mh0033/m300883/Tel_MMLE/docs/source/plots/monthly/month_extrc_diagram_patch_dual_time.png"
)
# %%
month_extrc_dig.bar_dual(Regress_slope,xlim = (-10,18))
plt.savefig(
    "/work/mh0033/m300883/Tel_MMLE/docs/source/plots/monthly/month_extrc_diagram_dual_tsurf.png"
)
# %%
