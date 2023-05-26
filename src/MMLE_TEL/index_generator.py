"""
To generate the index by projecting the data of all the years to the fixed spatial pattern.
"""

#%%
import xarray as xr
import numpy as np
import pandas as pd

import src.Teleconnection.vertical_eof as vertical_eof
import src.Teleconnection.tools as tools
import src.Teleconnection.rolling_eof as rolling_eof
import src.warming_stage.warming_stage as warming_stage
import src.Teleconnection.spatial_pattern as ssp
import src.MMLE_TEL.standardize as standardize


#%%
class decompose_troposphere:
    """
    A class to generate the eof and index of the whole troposphere
    """

    def __init__(
        self, model, vertical_eof, fixed_pattern="decade", standard="temporal_ens"
    ) -> None:
        self.vertical_eof = vertical_eof
        self.independence = self.vertical_eof == "ind"
        self.fixed_pattern = fixed_pattern
        self.model = model
        self.odir = "/work/mh0033/m300883/Tel_MMLE/data/" + self.model + "/"
        self.zg_path = self.odir + "zg_processed/"
        self.save_path = self.odir + "EOF_result/"
        self.standard = standard  # 'temporal', 'temporal_ens'

        # read data
        print("reading the gph data ...")
        self.data = read_data(self.zg_path)

        # decompose
        self.eof_result = self.decompose()
        self.std_eof_result = standard_index(self.eof_result, self.standard)

        # save
        self.save_result()

    def decompose(self):
        # deompose
        print("decomposing the data ...")
        eof_result = vertical_eof.vertical_eof(
            self.data,
            nmode=2,
            window=10,
            fixed_pattern=self.fixed_pattern,
            independent=self.independence,
        )
        return eof_result

    # save
    def save_result(self):
        print("saving the result ...")
        # save the result

        self.std_eof_result.to_netcdf(
            self.save_path
            + "troposphere_"
            + self.vertical_eof
            + "_"
            + self.fixed_pattern
            + "_"
            + self.standard
            + "_eof_result.nc"
        )


##########################################
class decompose_plev:
    """
    A class for decomposition of one single plev only.
    """

    def __init__(self, model, fixedPattern, plev=50000, standard="temporal") -> None:
        self.model = model
        self.plev = plev
        self.fixedPattern = fixedPattern  # warming or decade
        self.standard = standard

        self.odir = "/work/mh0033/m300883/Tel_MMLE/data/" + self.model + "/"
        self.zg_path = self.odir + "zg_processed/"
        self.ts_mean_path = self.odir + "ts_processed/ens_fld_year_mean.nc"
        self.save_path = self.odir + "EOF_result/"

        # read gph data
        print("reading the gph data ...")
        self.data = read_data(self.zg_path, plev=self.plev)

        # read ts_mean data if needed
        self.ts_mean = None
        if self.fixedPattern == "warming":
            self.ts_mean = warming_stage.read_tsurf_fldmean(self.ts_mean_path)

        # deompose
        self.eof_result = self.decompose()

        # standardize
        self.std_eof_result = standard_index(self.eof_result, self.standard)

    def decompose(self):
        """
        decompose the data
        """
        print("decomposing ...")

        eof_result = rolling_eof.rolling_eof(
            self.data, fixed_pattern=self.fixedPattern, ts_mean=self.ts_mean
        )
        return eof_result

    def save_result(self):
        print("saving the result ...")
        # save the unstandardized result
        self.eof_result.to_netcdf(
            self.save_path
            + "plev_"
            + str(self.plev)
            + "_"
            + self.fixedPattern
            + "_none_eof_result.nc"
        )
        # save the standardized result
        self.std_eof_result.to_netcdf(
            self.save_path
            + "plev_"
            + str(self.plev)
            + "_"
            + self.fixedPattern
            + "_"
            + self.standard
            + "_eof_result.nc"
        )


def read_data(zg_path, plev=None):
    """
    read data quickly
    """
    gph_dir = zg_path
    # read MPI_onepct data
    try:
        zg_data = xr.open_dataset(gph_dir + "allens_zg.nc")
        zg_data = tools.split_ens(zg_data)
    except FileNotFoundError:
        zg_data = xr.open_mfdataset(
            gph_dir + "*.nc", combine="nested", concat_dim="ens", join="override"
        )
    try:
        zg_data = zg_data.var156
    except AttributeError:
        zg_data = zg_data.zg

    # time to datetime
    try:
        zg_data["time"] = zg_data.indexes["time"].to_datetimeindex()
    except AttributeError:
        zg_data["time"] = pd.to_datetime(zg_data.time)

    # demean
    print(" demean the ensemble mean...")
    zg_ens_mean = zg_data.mean(dim="ens")
    zg_demean = zg_data - zg_ens_mean

    # select one altitude
    if plev is not None:
        print(" select the specific plev...")
        zg_demean = zg_demean.sel(plev=plev)

    return zg_demean


def standard_index(eof_result, standard="temporal"):
    print(f"standardizing the index with {standard} ...")
    # standarize the index with the tmeporal mean and std
    eof_result = eof_result.copy()
    if standard == "temporal":
        eof_result["pc"] = (
            eof_result["pc"] - eof_result["pc"].mean(dim="time")
        ) / eof_result["pc"].std(dim="time")

    elif standard == "temproal_ens":
        eof_result["pc"] = (
            eof_result["pc"] - eof_result["pc"].mean(dim=("time", "ens"))
        ) / eof_result["pc"].std(dim=("time", "ens"))

    return eof_result
