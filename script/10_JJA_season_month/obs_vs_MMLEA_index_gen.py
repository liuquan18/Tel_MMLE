# %%
import xarray as xr
import sys
from scipy.signal import detrend

# %%
import src.Teleconnection.spatial_pattern as ssp
import src.MMLE_TEL.index_generator as index_generator

#%%
import importlib
importlib.reload(ssp)
importlib.reload(index_generator)


# %%
# read data
def read_EAR_JJA_data(time = slice("1940", "2022")):
    odir = "/work/mh0033/m300883/Tel_MMLE/data/ERA5/"
    data_JJA = []
    for month in ["Jun", "Jul", "Aug"]:
        print(f"    reading the gph data of {month} ...")
        zg_path = odir + "zg_" + month + "/ts_era5_1940_2022.nc"

        # read data of single month
        zg_data = xr.open_dataset(zg_path)
        zg_data = zg_data.var129
        zg_data = zg_data / 9.80665  # convert to geopotential height (unit m)

        # select 500hPa
        zg_data = zg_data.sel(plev=50000)

        # detrend over every grid
        print("     detrending the data ...")
        arr = zg_data.values
        arr = detrend(arr, axis=0, type="linear")
        zg_data = zg_data.copy(data=arr)

        data_JJA.append(zg_data)

    # combine JJA
    data_JJA = xr.concat(data_JJA, dim="time")
    data_JJA = data_JJA.sortby("time")
    data_JJA = data_JJA.sel(time=time)

    return data_JJA


def read_LE_month_data(model, time=slice("1940", "2022")):
    zg_data = []
    for month in ["Jun", "Jul", "Aug"]:
        odir = f"/work/mh0033/m300883/Tel_MMLE/data/{model}/zg_{month}/"
        month_data = index_generator.read_data(odir, plev=50000)
        zg_data.append(month_data)
    zg_data = xr.concat(zg_data, dim="time")
    zg_data = zg_data.sortby("time")
    zg_data = zg_data.sel(time=time)
    return zg_data


# %%
def decompose_SNAO(model, time = slice("1940", "2022")):
    print(f"******{model}********")

    # read data
    print("----reading data----")
    if model == "ERA5":
        data_JJA = read_EAR_JJA_data(time=time)
        dim = 'time' # along time
    else:
        data_JJA = read_LE_month_data(model, time=time)
        data_JJA = data_JJA.stack(com=["time", "ens"])
        dim = 'com' # along ensemble and time
    
    # decompose
    print("----decomposing----")
    eof_result = ssp.doeof(
        data_JJA, nmode=2, dim=dim, standard="pc_temporal_std"
    ) 

    # save result
    print("----saving result----")
    sdir =  "/work/mh0033/m300883/Tel_MMLE/data/ERA5/EOF_result/"
    eof_result.to_netcdf(
        sdir + f"plev_50000_{time.start}_{time.stop}_{model}_all.nc"
    )

# %%
models = ['ERA5','MPI_GE','CanESM2','CESM1_CAM5','MK36','GFDL_CM3']
#%%
# get the mindex from keyboard
mindex = int(sys.argv[1])

#%%
model = models[mindex]
print(f"model: {model}")
decompose_SNAO(model, time = slice("1940", "2022"))
# %%

for model in models:
    decompose_SNAO(model, time = slice("1940", "2022"))
# %%