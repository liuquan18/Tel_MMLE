#%%
import xarray as xr
import pandas as pd
import numpy as np

import src.extreme.period_pattern_extreme as extreme
import proplot as pplt

import src.MMLE_TEL.extrc_tsurf as extrc_tsurf
import src.plots.extrc_tsurf_scatter as extrc_tsurf_scatter

#%%
import importlib
importlib.reload(extrc_tsurf)
importlib.reload(extrc_tsurf_scatter)

# %%
class extr_counts_tsurf:
    def __init__(
        self,
        models: list, # the name of the models
        ens: list=None,
        hlayers: int=50000, # which layer to compare
        average: bool = None
    ):

        # models
        self.models = models    
        self.hlayers = hlayers
        self.average = average
        self.ens = ens

        # extreme count, and mean surface
        self.extr_counts, self.ts_means = self.extrc_tsurf_all()
        ziped = zip(self.extr_counts, self.ts_means)
        self.extr_ts_pairs = list(ziped)

    def read_data(self, pc_dir, ts_dir):
        """read data of single model"""
        index = xr.open_dataset(pc_dir)
        ts = xr.open_dataset(ts_dir)
        try:
            t_mean = ts.mean(dim = 'ens')
        except ValueError:
            t_mean = ts

        try:
            t_mean = t_mean.tsurf
        except AttributeError:
            t_mean = t_mean.ts
                
        index = index.pc
        t_ano = t_mean - t_mean.isel(time = 0)
        t_ano = t_ano.squeeze()

        return index, t_ano


    def extrc_tsurf_all(self):

        odir = "/work/mh0033/m300883/Tel_MMLE/data/"

        extr_counts = []
        ts_means = []

        for i,model in enumerate(self.models):
            pc_dir = odir+model + "/EOF_result/ind_first_pc.nc"
            ts_dir = odir+model + "/ts_processed/tsurf_mean.nc"

            index, ts = self.read_data(pc_dir, ts_dir)

            extr_count, ts_mean = extrc_tsurf.decadal_extrc_tsurf(index, ts,hlayers = self.hlayers)   
            if self.average:
                extr_count = extr_count/self.ens[i]             
            extr_counts.append(extr_count)
            ts_means.append(ts_mean)            


        return extr_counts, ts_means


    def plot_scatter(self):
        extrc_tsurf_scatter.extCount_tsurf_scatter(self.extr_ts_pairs,labels = self.models)





# %%
allm = extr_counts_tsurf(['MPI_GE','MPI_GE_onepct','CanESM2'],ens = [100,100,40],hlayers = 20000,average=True)
# %%
