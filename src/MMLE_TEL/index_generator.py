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
        self.odir = "/work/mh0033/m300883/Tel_MMLE/data/" + self.model + "/"
        self.zg_path = (
            self.odir + "zg_processed/"
        )
        self.save_path = (
            self.odir + "EOF_result/"
        )

        # read data
        print("reading the gph data ...")
        self.data = season_eof.read_data(self.zg_path)
        self.eof_result = self.decompose()
        self.std_eof_result = self.standard_index()

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
        if (
            self.fixed_pattern == "first"
            or self.fixed_pattern == "last"
            or self.fixed_pattern == "all"
        ):
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

        else:
            print("wrong fixed pattern")
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


# %%
