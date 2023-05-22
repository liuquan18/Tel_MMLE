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
import src.MMLE_TEL.standardize as standardize


#%%
class decompose_fixedPattern:
    """
    A class to generate the eof and index
    """

    def __init__(
        self, model, vertical_eof, fixed_pattern="warming", standard="all"
    ) -> None:
        self.vertical_eof = vertical_eof
        self.independence = self.vertical_eof == "ind"
        self.fixed_pattern = fixed_pattern
        self.model = model
        self.odir = "/work/mh0033/m300883/Tel_MMLE/data/" + self.model + "/"
        self.zg_path = self.odir + "zg_processed/"
        self.save_path = self.odir + "EOF_result/"
        self.standard = standard  # 'own','all','none'

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
        if self.standard == "own":
            self.eof_result = standardize.standard_by_own(self.eof_result)
        elif self.standard == "all":
            try:
                all_index = xr.open_dataset(
                    self.odir + "EOF_result/" + self.vertical_eof + "_all_eof_result.nc"
                )

            except FileNotFoundError:
                print("all index not found, generate it first")
                # rase error
                raise FileNotFoundError

            self.eof_result = standardize.standard_by_all(all_index, self.eof_result)
        elif self.standard == "none":
            pass
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
            + self.standard
            + "_eof_result.nc"
        )


##########################################
class decompose_mmle:
    """
    A class for mmle decomposition
    """

    def __init__(self, model, fixedPattern, plev=50000, standarize=True) -> None:
        self.model = model
        self.plev = plev
        self.fixedPattern = fixedPattern  # warming or decade

        self.odir = "/work/mh0033/m300883/Tel_MMLE/data/" + self.model + "/"
        self.zg_path = self.odir + "zg_processed/"
        self.ts_mean_path = self.odir + "ts_processed/ens_fld_year_mean.nc"
        self.save_path = self.odir + "EOF_result/"

        # read gph data
        print("reading the gph data ...")
        self.data = self.read_data()

        # read ts_mean data if needed
        self.ts_mean = None
        if self.fixedPattern == "warming":
            self.ts_mean = warming_stage.read_tsurf_fldmean(self.ts_mean_path)

        # deompose
        self.eof_result = self.decompose()

        # standardize
        if standarize:
            self.std_eof_result = self.standard_index()

    def read_data(self):
        """
        read data quickly
        """
        gph_dir = self.zg_path
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
        print(" select the specific plev...")
        zg_gph = zg_demean.sel(plev=self.plev)

        return zg_gph

    def decompose(self):
        """
        decompose the data
        """
        print("decomposing ...")

        eof_result = rolling_eof.rolling_eof(
            self.data, fixed_pattern=self.fixedPattern, ts_mean=self.ts_mean
        )
        return eof_result

    def standard_index(self):
        print("standardizing the index ...")
        # standarize the index with the tmeporal mean and std
        eof_result = self.eof_result.copy()
        eof_result["pc"] = (
            eof_result["pc"] - eof_result["pc"].mean(dim="time")
        ) / eof_result["pc"].std(dim="time")
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
            + "_own_eof_result.nc"
        )
