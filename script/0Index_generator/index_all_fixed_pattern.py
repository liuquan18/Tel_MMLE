"""
To generate the index by projecting the data of all the years to the fixed spatial pattern.
"""

#%%
import xarray as xr
import numpy as np
import src.Teleconnection.season_eof as season_eof
import src.Teleconnection.tools as tools
import src.extreme.period_pattern_extreme as extreme

#%%
class decompose_fixedPattern:
    """
    A class to generate the eof and index
    """

    def __init__(self, model, vertical_eof, fixed_pattern, standard_ens_time=True):
        self.vertical_eof = vertical_eof
        self.independence = self.vertical_eof == "ind"
        self.fixed_pattern = fixed_pattern
        self.standard_ens_time = standard_ens_time
        self.model = model

        self.zg_path = "/work/mh0033/m300883/Tel_MMLE/data/"+self.model+"/zg_processed/*.nc"
        self.data = self.read_data()
        self.eof, self.index, self.fra = self.decompose()

        # standardize the index with mean and std of ens and time (150 yrs)
        if standard_ens_time:
            self.index = self.standardize()

        self.save_path = (
                        '/work/mh0033/m300883/Tel_MMLE/data/'+self.model+'/EOF_result/'
        )


    def read_data(self):
        zg_data = xr.open_mfdataset(self.zg_path, combine = 'nested',concat_dim = 'ens')
        # demean
        zg_ens_mean = zg_data.mean(dim = 'ens')
        zg_demean = zg_data-zg_ens_mean

        # select trop
        zg_trop = zg_demean.sel(plev = slice(100000,20000))
        return zg_trop

    def decompose(self):
        print("decomposing ...")
        eof_all, index_all, fra_all = season_eof.season_eof(
            self.data,
            nmode=2,
            window=10,
            fixed_pattern=self.fixed_pattern,
            independent=self.independence,
            standard=False,
        )
        return eof_all, index_all, fra_all

    def standardize(self, dim=("time", "ens")):
        """
        standardardize with the mean and std of 'time' and 'ens'.
        """
        print("standardizing...")
        mean = self.index.mean(dim=dim)
        std = self.index.std(dim=dim)
        index = (self.index - mean) / std
        return index

    def save_result(self):
        print("saving...")
        self.eof.to_netcdf(
            self.save_path
            + self.vertical_eof
            + "_"
            + self.fixed_pattern
            + "_"
            + "eof.nc"
        )
        self.index.to_netcdf(
            self.save_path
            + self.vertical_eof
            + "_"
            + self.fixed_pattern
            + "_"
            + "pc.nc"
        )
        self.fra.to_netcdf(
            self.save_path
            + self.vertical_eof
            + "_"
            + self.fixed_pattern
            + "_"
            + "fra.nc"
        )


# %%
cesm1_cam5_zg = xr.open_mfdataset("/work/mh0033/m300883/Tel_MMLE/data/CESM1_CAM5/zg_processed/*.nc")
# demean
cesm1_cam5_zg_ens_mean = cesm1_cam5_zg.mean(dim = 'ens')
cesm1_cam5_zg = cesm1_cam5_zg - cesm1_cam5_zg_ens_mean




#%%
ind_all = decompose_fixedPattern('ind','all')
ind_all.save_result()

#%%
dep_all = decompose_fixedPattern('dep','all')
dep_all.save_result()
# %%

ind_first = decompose_fixedPattern('ind','first')
ind_first.save_result()

# %%
dep_first = decompose_fixedPattern('dep','first')
dep_first.save_result()
# %%
