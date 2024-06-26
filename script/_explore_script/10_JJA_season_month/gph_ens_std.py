# %%
import xarray as xr
import numpy as np
import os
import sys

import src.MMLE_TEL.index_generator as index_generator
import src.MMLE_TEL.gph_statistic as gph_statistic
import src.Teleconnection.tools as tools

#%%
import importlib
importlib.reload(gph_statistic)


# %%
class gph_stats:
    def __init__(self,model,standard = 'first10') -> None:
        self.model = model
        self.standard = standard

        # dir to save the result
        self.to_dir = f"/work/mh0033/m300883/Tel_MMLE/data/{model}/box_based/"

        if not os.path.exists(self.to_dir):
            os.makedirs(self.to_dir)

        # # box statistics
        # print("-----------------------")
        # # self.boxStats = self.box_varibility()

        # slope of the ensemble std
        print("-----------------------")
        self.slopeStd = self.slope_ens_std()

        # slope of the ensemble mean
        print("-----------------------")
        self.slopeMean = self.slope_ens_mean()

        # # slope of the extreme events occurence
        # print("-----------------------")
        # self.slopeExtrc = self.slope_gph_extrc(standard = self.standard)

        # # slope of the box covariance
        # print("-----------------------")
        # self.slopeCov = self.slope_box_cov(standard = self.standard)

    
    def box_varibility(self):
        print(f"calculating the variability of {self.model} ...")
        model_pos, model_neg = gph_statistic.box_variability(self.model)
        # mkdir if not exist
        # create directory if it does not exist
        try:
            model_pos.to_netcdf(self.to_dir + "pos_var.nc")
            model_neg.to_netcdf(self.to_dir + "neg_var.nc")
        except PermissionError:
            os.remove(self.to_dir + "pos_var.nc")
            os.remove(self.to_dir + "neg_var.nc")
            model_pos.to_netcdf(self.to_dir + "pos_var.nc")
            model_neg.to_netcdf(self.to_dir + "neg_var.nc")

    def ens_std(self):
        print("calculating the ensemble std ...")
        std = gph_statistic.ens_std(self.model)
        # mkdir if not exist
        # create directory if it does not exist
        try:
            std.to_netcdf(self.to_dir + "ens_std.nc")
        except PermissionError:
            os.remove(self.to_dir + "ens_std.nc")
            std.to_netcdf(self.to_dir + "ens_std.nc")
    
    def slope_ens_std(self):
        print(f"calculating the slope of std {self.model} ...")
        slope = gph_statistic.slope_ens_std(self.model)
        # mkdir if not exist
        # create directory if it does not exist
        try:
            slope.to_netcdf(self.to_dir + "slope_of_ens_std.nc")
        except PermissionError:
            os.remove(self.to_dir + "slope_of_ens_std.nc")
            slope.to_netcdf(self.to_dir + "slope_of_ens_std.nc")

    
    def slope_ens_mean(self,**kwargs):
        print("calculating the slope of ens mean ...")
        slope = gph_statistic.slope_ens_mean(self.model,**kwargs)
        # mkdir if not exist
        # create directory if it does not exist

        try:
            slope.to_netcdf(self.to_dir + "slope_of_ens_mean.nc")
        except PermissionError:
            os.remove(self.to_dir + "slope_of_ens_mean.nc")
            slope.to_netcdf(self.to_dir + "slope_of_ens_mean.nc")

    def slope_gph_extrc(self,**kwargs):
        print("counting the extreme events ...")
        slope = gph_statistic.slope_gph_extrc(self.model,**kwargs)

        try:
            slope.to_netcdf(self.to_dir + f"slope_of_gph_extrc.nc")
        except PermissionError:
            os.remove(self.to_dir + f"slope_of_gph_extrc.nc")
            slope.to_netcdf(self.to_dir + f"slope_of_gph_extrc.nc")


    def slope_box_cov(self,**kwargs):
        print("calculating the slope of box covariance ...")
        slope = gph_statistic.slope_box_cov(self.model,**kwargs)

        try:
            slope.to_netcdf(self.to_dir + f"slope_of_box_cov.nc")
        except PermissionError:
            os.remove(self.to_dir + f"slope_of_box_cov.nc")
            slope.to_netcdf(self.to_dir + f"slope_of_box_cov.nc")

# %%
models = ["MPI_GE_onepct", "MPI_GE", "CanESM2", "CESM1_CAM5", "MK36", "GFDL_CM3"]

#%%
# get the num from keyboard
num = int(sys.argv[1])
# for model in models:
model = models[num - 1]

#%%
print(f"**********{model}**********")
gph_stats(model,standard='all')
