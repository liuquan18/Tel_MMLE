#%%
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
# %%
import src.reanalysis.remove_force as remove_force
import src.MMLE_TEL.index_generator as index_generate
import glob
#%%
import importlib
importlib.reload(remove_force)
#%%
def read_gph_data(model):
    odir = f"/work/mh0033/m300883/Tel_MMLE/data/{model}/"
    data_JJA = []
    for month in ["Jun", "Jul", "Aug"]:
        print(f"reading the gph data of {month} ...")
        zg_path = odir + "zg_" + month + "/"
        if model == 'MPI_GE':
            data_month = index_generate.read_data(zg_path, plev=50000,remove_ens_mean=False)
        elif model == 'CR20_allens':
            file_names = sorted(glob.glob(zg_path + "*.nc"))
            data_month = xr.open_mfdataset(file_names, combine="nested", concat_dim="ens")
            data_month = data_month['HGT']
        data_JJA.append(data_month)
    data = xr.concat(data_JJA, dim="time").sortby("time")
    return data

#%%
def _ens_std(data):
    return data.std()

def ens_std(data,ax, label, region = None, group_size = 10):
    # select the region and time
    # data = data.sel(time = slice('1850','2015'))
    if region is not None:
        data = data.sel(lon = slice(region[0],region[1]), lat = slice(region[2],region[3]))

    weights = np.cos(np.deg2rad(data.lat))
    weights.name = "weights"
    data_weighted = data.weighted(weights)
    glm = data_weighted.mean(dim = ('lon','lat'))
    glm_ens_std = glm.resample(time = f'{str(group_size)}Y').apply(_ens_std)
    glm_ens_std = glm_ens_std.sel(time = slice('1850','2015'))
    return glm_ens_std
    

# %%
# read 20CR data
CR20 = read_gph_data("CR20_allens")
MPI_GE = read_gph_data("MPI_GE")
# %%
# remove the forced signal
CR20_inter = remove_force.detrend(CR20, method="quadratic_trend")
MPI_GE_inter = remove_force.detrend(MPI_GE, method="quadratic_trend")
#%%
# plot the ensemble std
def compare_ens_std(region = None, save = False,group_size = 10):
    fig,ax = plt.subplots()
    glm_ens_std_CR20 = ens_std(CR20_inter,ax,label = '20CR', region = region, group_size = group_size)
    glm_ens_std_MPI_GE = ens_std(MPI_GE_inter,ax,label = 'MPI-GE', region = region, group_size = group_size)
    ax.set_title('Ensemble standard deviation')
    glm_ens_std_CR20.plot(ax = ax, label = '20CR', color = 'blue')
    glm_ens_std_MPI_GE.plot(ax = ax, label = 'MPI-GE', color = 'red')
    ax.legend()
    if save:
        plt.savefig(f'/work/mh0033/m300883/Tel_MMLE/docs/source/plots/supplyment/MPI_GE_20CR_ens_std.png')
# %%
compare_ens_std(group_size = 20)
# %%
compare_ens_std(group_size = 20,region = [-30,10,60,40],save=True)


# %%
compare_ens_std(group_size = 20,region = [-30,10,60,40])

# %%
