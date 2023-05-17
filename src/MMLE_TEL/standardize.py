# package to standardize indexes
#%%
import numpy as np
import pandas as pd
import xarray as xr

#%%

def standard_by_own(eof_result):
    print(" standardizing the index by its own temporal mean and std ...")
    eof_result["pc"] = (
        eof_result["pc"] - eof_result["pc"].mean(dim="time")
    ) / eof_result["pc"].std(dim="time")
    return eof_result

def standard_by_all(all_index, eof_result):
    print(" standardizing the index by temporal and ensemble mean and std ...")

    eof_result["pc"] = (
        eof_result["pc"] - all_index["pc"].mean(dim=("time", "ens"))
    ) / all_index["pc"].std(dim=("time", "ens"))
    return eof_result

