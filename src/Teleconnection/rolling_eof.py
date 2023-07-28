#%%
import numpy as np
import pandas as pd
import xarray as xr
from tqdm.notebook import tqdm, trange

import datetime
import src.Teleconnection.spatial_pattern as ssp
import src.Teleconnection.tools as tools
import src.warming_stage.warming_stage as warming_stage
import datetime

#%%
def rolling_eof(xarr, **kwargs):
    """do eof analysis for the data of every ten years.
    Parameters
    ----------
    xarr : xarray.DataArray
        the data to be decomposed.
    nmode : int, optional
        the number of mode to be generated, by default 2
    window : int, optional
        the window size if fix_pattern = 'False' is adopted., by default 10
    fixed_pattern : str, optional
        the fixed_pattern mode, by default "decade", can be "all", "warming" or "decade"
    ts_mean : xarray.DataArray, optional
        the mean of the time series, by default None
    Returns
    -------
    xarray.Dataset
        the result of eof analysis.
    """
    # get the parameters
    nmode = kwargs.get("nmode", 2)
    window = kwargs.get("window", 10)
    fixed_pattern = kwargs.get("fixed_pattern", "decade")
    ts_mean = kwargs.get("ts_mean", None)


    # if do the all-all decompose
    if fixed_pattern == "all":  # a little different from the following two.
        print("     using the all pattern")
        eof_result = ssp.doeof(
            tools.stack_ens(xarr, withdim="time"),
            nmode=nmode,
            dim="com",
            standard=False,
        )

    elif fixed_pattern == "warming":
        print("     decompose the warming period")

        # get the periods where the glmt increases 0K and 4K
        ts_mean = ts_mean
        warming_periods = warming_stage.temp_period(ts_mean)
        warming_index = xr.IndexVariable("warming", ["0K", "2K", "4K"])

        eof_results = []
        for period in warming_periods:
            print("         decomposing the warming period of {}".format(period))
            # decomose the decade
            eof_result_single = decompose_single_decade(xarr, period)
            eof_results.append(eof_result_single)

        eof_result = xr.concat(eof_results, dim=warming_index)

    elif fixed_pattern == "decade":
        print("     decomposing everty ten years")
        eof_result = decompose_decade(xarr, window)

    return eof_result

def decompose_decade(xarr, window):
    """decompose the data every ten years."""
    # start time
    years = np.unique(xarr.time.dt.year).astype('str')
    time_s = years[::window]
    # end time
    time_e = years[window-1::window]

    # create slice for each decade
    decade_slice = [slice(s, e) for s, e in zip(time_s, time_e)]

    # a list for storing the subarrays
    eofs = []
    pcs = []
    fras = []

    for time in decade_slice:
        print(f"     decomposing the decade of {time.start.astype(datetime64[Y])} - {time.stop.astype(datetime64[Y])}")
        # slice the time

        eof_result_single = decompose_single_decade(xarr, time)
        eof = eof_result_single["eof"]
        pc = eof_result_single["pc"].copy()
        fra = eof_result_single["fra"]

        eofs.append(eof)
        pcs.append(pc)
        fras.append(fra)

    # concat the subarrays together, and make the decade as a new dim
    EOF = xr.concat(eofs, dim="decade")
    FRA = xr.concat(fras, dim="decade")
    PC = xr.concat(pcs, "time")

    # combine EOF, FRA, PC together as a dataset
    eof_result = xr.Dataset({"eof": EOF, "pc": PC, "fra": FRA})
    return eof_result


def decompose_single_decade(xarr, timeslice, nmode=2):
    """decompose a single decade."""
    field = xarr.sel(time=timeslice)
    field = field.stack(com=("ens", "time"))

    eof_result = ssp.doeof(field, nmode=nmode, dim="com")

    return eof_result


