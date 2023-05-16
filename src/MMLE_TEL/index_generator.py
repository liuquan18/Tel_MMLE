"""
To generate the index by projecting the data of all the years to the fixed spatial pattern.
"""

#%%
import xarray as xr
import numpy as np
import pandas as pd

import src.Teleconnection.season_eof as season_eof
import src.Teleconnection.tools as tools
import src.Teleconnection.rolling_eof as rolling_eof
import src.warming_stage.warming_stage as warming_stage
import src.Teleconnection.spatial_pattern as ssp


#%%
class decompose_fixedPattern:
    """
    A class to generate the eof and index
    """

    def __init__(self, model, vertical_eof, fixed_pattern = 'warming'):
        self.vertical_eof = vertical_eof
        self.independence = self.vertical_eof == "ind"
        self.fixed_pattern = fixed_pattern
        self.model = model
        self.odir = "/work/mh0033/m300883/Tel_MMLE/data/" + self.model + "/"
        self.zg_path = self.odir + "zg_processed/"
        self.save_path = self.odir + "EOF_result/"

        # read data
        print("reading the gph data ...")
        self.data = season_eof.read_data(self.zg_path)
        self.eof_result = self.decompose()
        self.std_eof_result = self.standard_index()
        self.save_result()

    def decompose(self):
        # deompose
        print("decomposing the data ...")
        eof_result = season_eof.season_eof(
            self.data,
            nmode=2,
            window=10,
            fixed_pattern=self.fixed_pattern,
            independent=self.independence,
            standard=False,  # standardization is done seperately
        )
        return eof_result

    def standard_index(self):
        print("standardizing the index ...")
        # if single pattern, standardize the index with its own mean and std (temporal and ens)
        if self.fixed_pattern == "first" or self.fixed_pattern == "last":
            self.eof_result["pc"] = (
                self.eof_result["pc"] - self.eof_result["pc"].mean(dim=("time", "ens"))
            ) / self.eof_result["pc"].std(dim=("time", "ens"))

        # if changing pattern, standardize the index with the mean and std of all the index
        elif self.fixed_pattern == "decade" or self.fixed_pattern == "False":
            try:
                all_index = xr.open_dataset(
                    self.odir + "EOF_result/" + self.vertical_eof + "_all_eof_result.nc"
                )

            except FileNotFoundError:
                print("all index not found, generate it first")

            self.eof_result["pc"] = (
                self.eof_result["pc"] - all_index["pc"].mean(dim=("time", "ens"))
            ) / all_index["pc"].std(dim=("time", "ens"))

        elif self.fixed_pattern == "all":
            print("     no standarization for all pattern")
        return self.eof_result

    # save
    def save_result(self):
        print("saving the result ...")
        # save the result

        self.std_eof_result.to_netcdf(
            self.save_path
            + self.vertical_eof
            + "_"
            + self.fixed_pattern
            + "_"
            + "eof_result.nc"
        )


##########################################
class decompose_mmle:
    """
    A class for mmle decomposition
    """

    def __init__(self, model, fixedPattern, gph=50000) -> None:
        self.model = model
        self.gph = gph
        self.fixedPattern = fixedPattern  # warming or decade

        self.odir = "/work/mh0033/m300883/Tel_MMLE/data/" + self.model + "/"
        self.zg_path = self.odir + "zg_processed/"
        self.ts_mean_path = self.odir + "ts_processed/ens_fld_year_mean.nc"
        self.save_path = self.odir + "EOF_result/"

        # read data
        print("reading the gph data ...")
        self.data = self.read_data()

        # warming stages
        if self.fixedPattern == "warming":
            self.ts_mean = warming_stage.read_tsurf_fldmean(self.ts_mean_path)
            self.warming_periods = warming_stage.temp_period(self.ts_mean)

            # decompose
            self.all_eof, self.eof_0K, self.eof_4K = self.decompose_warming_eof()
        
        elif self.fixedPattern == "decade":
            self.decade_eof = self.decompose_decade_eof()

    def read_data(self):
        """
        read data quickly
        """
        gph_dir = self.zg_path
        # read MPI_onepct data
        try:
            zg_data = xr.open_dataset(gph_dir + "allens.nc")
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
        print(" select the specific gph...")
        zg_gph = zg_demean.sel(plev=self.gph)
        # standardize seperately with the temporal mean and std
        print(" standardize each altitudes seperately...")
        zg_gph = (zg_gph - zg_gph.mean(dim="time")) / zg_gph.std(dim="time")
        return zg_gph

    def decompose_allPattern(self):
        """
        decompose with the all pattern
        """
        all_eof = rolling_eof.rolling_eof(
            self.data, nmode=2, window=6, fixed_pattern="all"
        )
        return all_eof

    def decompose_warming_period_pattern(self, warming_stage):
        """
        decompose the data in the 10 years around the 0K or 4K warming stage
        """
        if warming_stage == "0K":
            warming_period = self.warming_periods[0]
        elif warming_stage == "4K":
            warming_period = self.warming_periods[-1]

        print(f"    decomposing the warming stage {warming_stage} ...")
        field = self.data.sel(time=warming_period)
        field = field.stack(com=("ens", "time"))
        eof_result = ssp.doeof(field, nmode=2, dim="com")
        eof_result = eof_result[["eof", "fra"]]
        PC = ssp.project_field(field, eof_result.eof, dim="com")
        eof_result["pc"] = PC
        eof_result.attrs["warming_stage"] = str(warming_period)
        return eof_result

    def decompose_warming_eof(self):
        """
        decompose and then standardize the index
        """
        print("decomposing the eofs ...")

        all_eof = self.decompose_allPattern()
        eof_0K = self.decompose_warming_period_pattern("0K")
        eof_4K = self.decompose_warming_period_pattern("4K")
        # standardize the index with the mean and std of the all index
        print("standardize the index ...")
        eof_0K["pc"] = (
            eof_0K["pc"] - all_eof["pc"].mean(dim=("time", "ens"))
        ) / all_eof["pc"].std(dim=("time", "ens"))
        eof_4K["pc"] = (
            eof_4K["pc"] - all_eof["pc"].mean(dim=("time", "ens"))
        ) / all_eof["pc"].std(dim=("time", "ens"))
        return all_eof, eof_0K, eof_4K


    def decompose_decade_eof(self):
        """
        decompose the historical and ssp8.5 every ten years.
        """
        print("decomposing every ten years ...")

        try:
            all_eof = xr.open_dataset(self.save_path + "gph_" + str(self.gph) + "_all_" + "eof_result.nc")
        except FileNotFoundError:
            all_eof = self.decompose_allPattern()

        decade_eof = rolling_eof.rolling_eof(
            self.data, nmode=2, window=10, fixed_pattern="decade",standard=False
        )
        print("standardize the index ...")
        decade_eof["pc"] = (
            decade_eof["pc"] - all_eof["pc"].mean(dim=("time", "ens"))
        ) / all_eof["pc"].std(dim=("time", "ens"))

        return decade_eof


    def save_result(self):
        print("saving the result ...")
        # save the result

        self.all_eof.to_netcdf(
            self.save_path + "gph_" + str(self.gph) + "_all_" + "eof_result.nc"
        )

        if self.fixedPattern == 'warming':
            self.eof_0K.to_netcdf(
                self.save_path + "gph_" + str(self.gph) + "_0K_" + "eof_result.nc"
            )
            self.eof_4K.to_netcdf(
                self.save_path + "gph_" + str(self.gph) + "_4K_" + "eof_result.nc"
            )
        elif self.fixedPattern == 'decade':
            self.decade_eof.to_netcdf(
                self.save_path + "gph_" + str(self.gph) + "_decade_" + "eof_result.nc"
            )
