"""
To generate the index by projecting the data of all the years to the fixed spatial pattern.
"""

# %%
import xarray as xr
import numpy as np
import pandas as pd
import random
import os

import src.Teleconnection.vertical_eof as vertical_eof
import src.Teleconnection.tools as tools
import src.Teleconnection.rolling_eof as rolling_eof
import src.warming_stage.warming_stage as warming_stage
import warnings
import glob


# %%
class decompose_troposphere:
    """
    A class to generate the eof and index of the whole troposphere
    """

    def __init__(
        self,
        model,
        vertical_eof="ind",
        fixedPattern="decade",
        standard="temporal_ens",
        season="DJFM",
        all_years=False,
    ) -> None:
        self.vertical_eof = vertical_eof
        self.independence = self.vertical_eof == "ind"
        self.fixed_pattern = fixedPattern
        self.model = model
        self.odir = "/work/mh0033/m300883/Tel_MMLE/data/" + self.model + "/"
        self.save_path = self.odir + "EOF_result/"
        self.standard = standard  # 'temporal', 'temporal_ens'
        self.season = season
        self.zg_path = self.odir + "zg_" + self.season + "/"

        # read data
        print(f"reading the gph data of {self.season} ...")
        # read gph data
        data_JJA = []
        for month in ["Jun", "Jul", "Aug"]:
            print(f"reading the gph data of {month} ...")
            zg_path = self.odir + "zg_" + month + "/"
            data_JJA.append(read_data(zg_path))
        data = xr.concat(data_JJA, dim="time").sortby("time")
        
        if all_years:
            self.data = data
        else:
            data_first = self.select_year(data, 0, 10)
            data_last = self.select_year(data, -20,-10)  # since the last 10 years is not complete (no data in 2100 in MPI_GE, no data in 2000 in MPI_GE_onepct)
            # also to keep the time range the same as the decompose_plev
            self.data = xr.concat([data_first, data_last], dim="time")

        # decompose
        self.eof_result = self.decompose()
        self.std_eof_result = standard_index(self.eof_result, self.standard)

        # save
        self.save_result()

    def decompose(self):
        # deompose
        print("decomposing the data over the whole troposhpere ...")
        eof_result = vertical_eof.vertical_eof(
            self.data,
            nmode=2,
            window=10,
            fixed_pattern=self.fixed_pattern,
            independent=self.independence,
        )
        return eof_result
    

    def select_year(self, data, year_index_start, year_index_end):
        """
        select the data of the specific years
        """
        years = np.unique(data["time.year"])
        years = sorted(years)[year_index_start:year_index_end]
        res = data.sel(time=data["time.year"].isin(years))
        return res
    # save
    def save_result(self):
        print("saving the result ...")
        # save the result
        self.eof_result.to_netcdf(
            self.save_path
            + "troposphere_"
            + self.vertical_eof
            + "_"
            + self.fixed_pattern
            + "_"
            + self.season
            + "_none_eof_result.nc"
        )

        self.std_eof_result.to_netcdf(
            self.save_path
            + "troposphere_"
            + self.vertical_eof
            + "_"
            + self.fixed_pattern
            + "_"
            + self.standard
            + "_"
            + self.season
            + "_eof_result.nc"
        )


# %%
##########################################
class decompose_plev:
    """
    A class for decomposition of one single plev only.
    """

    def __init__(
        self, model, fixedPattern, plev=50000, standard="first", season="DJFM"
    ) -> None:
        self.model = model
        self.plev = plev
        self.fixedPattern = fixedPattern  # warming or decade
        self.standard = standard
        self.season = season

        self.odir = "/work/mh0033/m300883/Tel_MMLE/data/" + self.model + "/"
        self.ts_mean_path = self.odir + "ts_processed/ens_fld_year_mean.nc"
        self.save_path = self.odir + "EOF_result/"
        self.zg_path = self.odir + "zg_" + self.season + "/"

        # read gph data
        print(f"reading the gph data of {self.season} ...")
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
        nondir = (
            self.save_path
            + "plev_"
            + str(self.plev)
            + "_"
            + self.fixedPattern
            + "_"
            + self.season
            + "_none_eof_result.nc"
        )
        stddir = (
            self.save_path
            + "plev_"
            + str(self.plev)
            + "_"
            + self.fixedPattern
            + "_"
            + self.standard
            + "_"
            + self.season
            + "_eof_result.nc"
        )
        try:
            self.eof_result.to_netcdf(nondir)
        except PermissionError:
            os.remove(nondir)
            self.eof_result.to_netcdf(nondir)
        try:
            self.std_eof_result.to_netcdf(stddir)
        except PermissionError:
            os.remove(stddir)
            self.std_eof_result.to_netcdf(stddir)


# %%
##########################################
class decompose_plev_random_ens:
    def __init__(
        self,
        fixedPattern,
        ens_size,
        base_model="MPI_GE_onepct",
        plev=50000,
        standard="first",
    ) -> None:
        self.model = base_model + "_random"
        self.plev = plev
        self.fixedPattern = fixedPattern  # warming or decade
        self.standard = standard
        self.ens_size = ens_size
        self.season = 'JJA'

        self.odir = "/work/mh0033/m300883/Tel_MMLE/data/" + base_model + "/"

        self.ts_mean_path = self.odir + "ts_processed/ens_fld_year_mean.nc"
        self.save_path = "/work/mh0033/m300883/Tel_MMLE/data/" + self.model + "/EOF_result/"

        # read gph data
        data_JJA = []
        for month in ["Jun", "Jul", "Aug"]:
            print(f"reading the gph data of {month} ...")
            zg_path = self.odir + "zg_" + month + "/"
            data_JJA.append(
                read_data(zg_path, plev=self.plev, remove_ens_mean=False)
            )  # not remove the ensemble mean here yet.
        self.data = xr.concat(data_JJA, dim="time").sortby("time")

        # randomly select ens_size members
        random.seed(1)
        self.data = self.data.isel(ens=random.sample(range(0, 100), self.ens_size))
        # remove the ensemble mean
        print("removing the ensemble mean ...")
        self.data = self.data - self.data.mean(dim="ens")

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

        eof_result = rolling_eof.rolling_eof(self.data, fixed_pattern=self.fixedPattern)
        return eof_result

    def save_result(self):
        print("saving the result ...")
        none_standard_name = (
            self.save_path
            + "plev_"
            + str(self.plev)
            + "_"
            + self.fixedPattern
            + "_"
            + self.season
            + "_"
            + str(self.ens_size)
            + "_none_eof_result.nc"
        )

        standard_name = (
                self.save_path
                + "plev_"
                + str(self.plev)
                + "_"
                + self.fixedPattern
                + "_"
                + self.season
                + "_"
                + self.standard
                + "_"
                + str(self.ens_size)
                + "_eof_result.nc"
        )

        try:
            # save the unstandardized result
            self.eof_result.to_netcdf(none_standard_name)
            # save the standardized result
            self.std_eof_result.to_netcdf(standard_name)
        except PermissionError:
            os.remove(none_standard_name)            
            os.remove(standard_name)
            # save the unstandardized result
            self.eof_result.to_netcdf(none_standard_name)
            # save the standardized result
            self.std_eof_result.to_netcdf(standard_name)


# %%
##########################################
class decompose_plev_JJA:
    """
    A class for decomposition of one single plev only.
    """

    def __init__(
        self,
        model,
        fixedPattern="decade",
        plev=50000,
        standard="first",
    ) -> None:
        self.model = model
        self.plev = plev
        self.fixedPattern = fixedPattern  # warming or decade
        self.standard = standard
        self.season = "JJA"

        self.odir = "/work/mh0033/m300883/Tel_MMLE/data/" + self.model + "/"
        self.ts_mean_path = self.odir + "ts_processed/ens_fld_year_mean.nc"
        self.save_path = self.odir + "EOF_result/"

        # read gph data
        data_JJA = []
        for month in ["Jun", "Jul", "Aug"]:
            print(f"reading the gph data of {month} ...")
            zg_path = self.odir + "zg_" + month + "/"
            data_JJA.append(read_data(zg_path, plev=self.plev))
        self.data = xr.concat(data_JJA, dim="time").sortby("time")

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
        nondir = (
            self.save_path
            + "plev_"
            + str(self.plev)
            + "_"
            + self.fixedPattern
            + "_"
            + self.season
            + "_none_eof_result.nc"
        )
        stddir = (
            self.save_path
            + "plev_"
            + str(self.plev)
            + "_"
            + self.fixedPattern
            + "_"
            + self.standard
            + "_"
            + self.season
            + "_eof_result.nc"
        )
        try:
            self.eof_result.to_netcdf(nondir)
        except PermissionError:
            os.remove(nondir)
            self.eof_result.to_netcdf(nondir)
        try:
            self.std_eof_result.to_netcdf(stddir)
        except PermissionError:
            os.remove(stddir)
            self.std_eof_result.to_netcdf(stddir)

# %%
##########################################

class projected_index:



    def __init__(
        self,
        model,
        plev=50000,
        standard="first",
    ) -> None:
        self.model = model
        self.plev = plev
        self.standard = standard
        self.season = "JJA"

        self.odir = "/work/mh0033/m300883/Tel_MMLE/data/" + self.model + "/"
        self.eofdir = "/work/mh0033/m300883/Tel_MMLE/data/" + self.model + "/EOF_result/"
        self.save_path = self.odir + "EOF_result/"

        # read gph data
        data_JJA = []
        for month in ["Jun", "Jul", "Aug"]:
            print(f"reading the gph data of {month} ...")
            zg_path = self.odir + "zg_" + month + "/"
            data_JJA.append(read_data(zg_path, plev=self.plev))
        self.data = xr.concat(data_JJA, dim="time").sortby("time")

        # read eof data
        eof = xr.open_dataset()


        # deompose
        self.eof_result = self.decompose()

        # standardize
        self.std_eof_result = standard_index(self.eof_result, self.standard)



# %%
##########################################


def read_data(
    zg_path,
    plev=None,
    remove_ens_mean=True,
    var_name = "zg",
):
    """
    read data quickly
    """
    gph_dir = zg_path
    # read MPI_onepct data
    # fix the order of ensemble members
    print("reading the gph data of all ensemble members...")
    all_ens_lists = sorted(
        glob.glob(gph_dir + "*.nc")
    )  # to make sure that the order of ensemble members is fixed
    zg_data = xr.open_mfdataset(
        all_ens_lists, combine="nested", concat_dim="ens", join="override",
    )  # consider chunks={}, # kz the file size is small (< 3G). 

    zg_data["ens"] = np.arange(zg_data.ens.size)
    try:
        zg_data = zg_data.var156
    except AttributeError:
        try:
            zg_data = zg_data.zg
    
        except AttributeError:
            zg_data = zg_data[var_name]

    # time to datetime
    try:
        zg_data["time"] = zg_data.indexes["time"].to_datetimeindex()
    except AttributeError:
        zg_data["time"] = pd.to_datetime(zg_data.time)

    # demean
    if remove_ens_mean:
        print(" demean the ensemble mean...")
        zg_ens_mean = zg_data.mean(dim="ens")
        zg_demean = zg_data - zg_ens_mean
    else:
        zg_demean = zg_data

    # select one altitude
    try:
        if plev is not None:
            print(" select the specific plev...")
            zg_plev = zg_demean.sel(plev=plev)
        else:
            # select the 1000hPa - 200hPa
            print(" select the 1000hPa - 200hPa...")
            zg_plev = zg_demean.sel(plev=slice(100000, 20000))
            if zg_plev.plev.size == 0:
                zg_plev = zg_demean.sel(plev=slice(20000, 100000))
    except KeyError:
        zg_plev = zg_demean # for the data only with one plev
    return zg_plev


def standard_index(eof_result, standard):
    # standarize the index with the tmeporal mean and std
    eof_result = eof_result.copy()
    if standard == "first":
        print(" standardizing the index with the first 10 years ...")
        years = np.unique(eof_result["time.year"])
        years = sorted(years)[:10]
        ref = eof_result["pc"].sel(time=eof_result["time.year"].isin(years))
        pc_std = (eof_result["pc"] - ref.mean(dim=("time", "ens"))) / ref.std(
            dim=("time", "ens")
        )
        eof_result["pc"] = pc_std

    # standarize the index with the temporal and ensemble mean and std
    elif standard == "temporal_ens":
        print(" standardizing the index with temporal and ensemble mean and std ...")
        eof_result = eof_result.copy()
        ref = eof_result["pc"]
        pc_std = (eof_result["pc"] - ref.mean(dim=("time", "ens"))) / ref.std(
            dim=("time", "ens")
        )
        eof_result["pc"] = pc_std
    return eof_result
