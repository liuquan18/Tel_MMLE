# %%
import numpy as np
import importlib
import matplotlib.pyplot as plt
import src.plots.extreme_plot as extrc_tsurf
import src.MMLE_TEL.index_stats as index_stats
import xarray as xr
import os
import sys
import src.composite as composite
#%%
import importlib

importlib.reload(index_stats)
importlib.reload(extrc_tsurf)
importlib.reload(composite)

#%%
def read_var_months(model,var_name='ts',remove_ensmean = True, **kwargs):
    odir = f'/work/mh0033/m300883/Tel_MMLE/data/{model}/'
    Jun_fl = odir + f'{var_name}_Jun/'
    Jul_fl = odir + f'{var_name}_Jul/'
    Aug_fl = odir + f'{var_name}_Aug/'

    Jun_f = index_stats.read_var_data(Jun_fl,remove_ensmean=remove_ensmean)
    Jul_f = index_stats.read_var_data(Jul_fl,remove_ensmean=remove_ensmean)
    Aug_f = index_stats.read_var_data(Aug_fl,remove_ensmean=remove_ensmean)

    JJA_f = xr.concat([Jun_f, Jul_f, Aug_f], dim='time')
    JJA_f = JJA_f.sortby('time')
    if var_name == 'ts':
        try:
            JJA_f = JJA_f['tsurf']
        except KeyError:
            JJA_f = JJA_f['ts']
    elif var_name == 'pr':
        try:
            JJA_f = JJA_f['pr']
        except KeyError:
            JJA_f = JJA_f['precip']
    elif var_name == 'zg':
        try:
            JJA_f = JJA_f.var156
        except AttributeError:
            JJA_f = JJA_f.zg
    elif var_name == 'u10':
        JJA_f = JJA_f.u10
    elif var_name == 'v10':
        JJA_f = JJA_f.v10
    elif var_name == 'bDays':
        JJA_f = JJA_f['IB index']
        JJA_f = JJA_f.drop('plev')

    elif var_name == 'u':
        JJA_f = JJA_f.var131


    return JJA_f
#%%
def same_count(model,fixed_pattern = 'decade_mpi'):
    """read extreme counts"""
    print("reading the extreme counts ... ")
    odir = "/work/mh0033/m300883/Tel_MMLE/data/" + model + "/extreme_count/"
    filename = f"plev_50000_{fixed_pattern}_first_JJA_extre_counts.nc"
    ds = xr.open_dataset(odir + filename).pc
    count = ds.isel(time = 0).sel(confidence = 'true')
    count = count.drop(['plev','confidence','time'])
    return count
#%%
def composite(model,var_name='ts',reduction = 'mean',**kwargs):
    remove_ensmean = kwargs.get('remove_ensmean',True)
    var_data = read_var_months(model,var_name=var_name,remove_ensmean=remove_ensmean)
    if var_name == 'zg':
        var_data = var_data.sel(plev=50000).drop('plev')
    count = same_count(model) if reduction == 'mean_same_number' else 'all'
    index_stats.composite_analysis(
        model            = model,
        index_season     = 'JJA',
        var_season     = 'JJA',
        fixed_pattern    = 'decade_mpi',
        var_data         = var_data,
        var_name         = var_name,
        reduction        = reduction,
        count            = count,
        )
    
#%%
num = int(sys.argv[1])
t1 = int(sys.argv[2])
t2 = int(sys.argv[3])

#%%
def mean_all(num,rank):
    models = ['MPI_GE_onepct','MK36','GFDL_CM3','CanESM2','CESM1_CAM5','MPI_GE']
    vars = ['zg']

    print("===========================================")
    print(f"node_num:{num} is doing {models[num-1]}")
    composite(models[num-1])

    # for reduction = 'mean'
    model = models[num-1]
    var_name = vars[rank]
    print("===========================================")

    print(f"model {model} var {var_name} is doing")
    composite(model,var_name=var_name)


#%%

def mean_same_number(num):
    models = ['MPI_GE_onepct','MK36','GFDL_CM3','CanESM2','CESM1_CAM5','MPI_GE']
    composite(models[num-1],reduction= 'mean_same_number')

#%%
# main run mean_all
if __name__ == '__main__':

    # === mpi4py ===
    try:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        npro = comm.Get_size()
    except:
        print('::: Warning: Proceeding without mpi4py! :::')
        rank = 0
        npro = 1

    mean_same_number(num)

# # %%
# composite('MPI_GE_onepct_30','bDays',remove_ensmean=False)
# # %%
# composite('MPI_GE_onepct_30','bDays_ano',remove_ensmean=False)

# # %%
# composite('GFDL_CM3',reduction = 'mean_same_number')
# # %%
