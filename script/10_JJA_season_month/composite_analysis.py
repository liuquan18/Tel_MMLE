# %%
import numpy as np
import importlib
import matplotlib.pyplot as plt
import src.plots.extreme_plot as extrc_tsurf
import src.MMLE_TEL.index_stats as index_stats
import xarray as xr
import os
import sys
#%%
import importlib

importlib.reload(index_stats)
importlib.reload(extrc_tsurf)

# %%
import multiprocessing as mp
from multiprocessing import Pool

#%%
def read_var_months(model,var_name='ts'):
    odir = f'/work/mh0033/m300883/Tel_MMLE/data/{model}/'
    Jun_fl = odir + f'{var_name}_Jun/'
    Jul_fl = odir + f'{var_name}_Jul/'
    Aug_fl = odir + f'{var_name}_Aug/'

    Jun_f = index_stats.read_var_data(Jun_fl)
    Jul_f = index_stats.read_var_data(Jul_fl)
    Aug_f = index_stats.read_var_data(Aug_fl)

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

    return JJA_f

def composite(model,var_name='ts',reduction = 'mean',count = None):
    var_data = read_var_months(model,var_name=var_name)
    index_stats.composite_analysis(
        model            = model,
        index_season     = 'JJA',
        var_season     = 'JJA',
        fixed_pattern    = 'decade_mpi',
        var_data         = var_data,
        var_name         = var_name,
        reduction        = reduction,
        count            = count
        )
    
#%%
num = int(sys.argv[1])
t1 = int(sys.argv[2])
t2 = int(sys.argv[3])

#%%
def mean_all(num,rank):
    models = ['MPI_GE_onepct','MK36','GFDL_CM3','CanESM2','CESM1_CAM5','MPI_GE']
    vars = ['ts','pr']

    print("===========================================")
    print(f"node_num:{num} is doing {models[num-1]}")
    composite(models[num-1])

    # for reduction = 'mean'
    model = models[num-1]
    var_name = vars[rank]
    print("===========================================")

    print(f"model {model} var {var_name} is doing")
    composite(model,var_name=var_name)

def mean_same_number():
    # for reduction = 'mean_same_number'
    composite('CESM1_CAM5',reduction= 'mean_same_number',count = 80)
    composite('MK36',reduction= 'mean_same_number',count = 50)
    composite('CanESM2',reduction= 'mean_same_number',count = 80)
    composite('GFDL_CM3',reduction= 'mean_same_number',count = 30)
    composite('MPI_GE_onepct',reduction= 'mean_same_number',count = 200)
    composite('MPI_GE',reduction= 'mean_same_number',count = 200)

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

    mean_all(num,rank)

#%%
mean_same_number()
# %%
