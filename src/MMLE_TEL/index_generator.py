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

        self.zg_path = (
            "/work/mh0033/m300883/Tel_MMLE/data/" + self.model + "/zg_processed/"
        )

        # read data
        print("reading the gph data ...")
        self.data = season_eof.read_data(self.zg_path)

        # deompose
        print("decomposing the data ...")
        self.eof_result = season_eof.season_eof(
            self.data,
            nmode=2,
            window=10,
            fixed_pattern=self.fixed_pattern,
            independent=self.independence,
            standard=False,
        )

        # save the result
        self.save_path = (
            "/work/mh0033/m300883/Tel_MMLE/data/" + self.model + "/EOF_result/"
        )
        self.eof_result.to_netcdf(
            self.save_path
            + self.vertical_eof
            + "_"
            + self.fixed_pattern
            + "_"
            + "eof_result.nc"
        )

 

# %%
