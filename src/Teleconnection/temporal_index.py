"""
This is the source code for temporal index generator
"""

import numpy as np
import xarray as xr

import src.Teleconnection.season_eof as season_eof
import src.Teleconnection.spatial_pattern as ssp


def index_diff_pattern(xarr, independent=True, standard=True):
    """
    projeting the whole time series onto the three different patterns to get
    the corresponding index.
    Three patterns:
        - all: all the ensembles and years are fed to the EOF.
        - first: the ensembles of the first 10 years are fed to EOF.
        - last: the ensembles of the last 10 years are fed to EOF.
    The above three patterns can be derived with all-level-independent or
    all-level-dependent, the latter is a common pattern for all levels.
    **Argument**:
        *xarr* the data to be projected and to get the spatial patterns
        *independent* how to get the spatil patterns.
        *standard* whether to do the standarization or not. here all the
                   three indexes are standardized with the temporal mean
                   and std of index from all-pattern.
    **Return**:
        Three indexes projected from the same time series (all the years) but onto
        three different spatial patterns.
    """
    _, all_all, _ = season_eof.season_eof(
        xarr, nmode=2, window=10, fixed_pattern="all", independent=independent
    )  # "method" doesn't matter since the pc is
    # calculated independently.
    _, all_first, _ = season_eof.season_eof(
        xarr, nmode=2, window=10, fixed_pattern="first", independent=independent
    )
    _, all_last, _ = season_eof.season_eof(
        xarr, nmode=2, window=10, fixed_pattern="last", independent=independent
    )

    # using the mean and std of all-all to standardize all the series.
    if standard:
        all_mean = all_all.mean(dim="time")
        all_std = all_all.std(dim="time")

        all_all = (all_all - all_mean) / all_std
        all_first = (all_first - all_mean) / all_std
        all_last = (all_last - all_mean) / all_std

    return all_all, all_first, all_last


def index_changing_pattern(xarr, independent=True, standard=True):
    """
    generate the index from the changing pattern.
    projecting the data series (10 years and 100 ensembles) onto the corresponding
    spatial patterns, which varies as time goes by. such spatial patterns are generated
    by so called 'rolling_eof'.
    """
    EOF, index, FRA = season_eof.season_eof(
        xarr,
        nmode=2,
        method="rolling_eof",
        window=10,
        fixed_pattern=False,
        return_full_eof=True,
        independent=independent,
        standard=True,
    )

    if standard:
        index_mean = index.mean(dim="time")
        index_std = index.std(dim="time")
        index = (index - index_mean) / index_std

    return EOF, index, FRA


def period_index(all_indexes, period):
    """
    select the first10 or last10 index from the whole time series of teleconnction index.
    the index name can be seen from the table below:
    |spatial pattern|   all   |    first    |    last   |
    |temporal period|
    |---------------|---------|-------------| ----------|
    |first 10       |first10-all| first10-first | first10-last|
    |last 10        |last10-all | last10-first  | last10-last |
    **Arguments**
        *all_index* the three index of all years to the three patterns.
                    should be order in [all_all, all_first, all_last]
        *period* the first 10 or last 10 years of index
    **Return**
        the three index for first or last 10 years. ordered in
        [_all,_first,_last]
    """
    if period == "first10":
        ten_indexes = [all_index.isel(time=slice(0, 10)) for all_index in all_indexes]
    if period == "last10":
        ten_indexes = [
            all_index.isel(time=slice(-10, None)) for all_index in all_indexes
        ]
    return ten_indexes

def main():
    all_all_ind = xr.open_dataset(
        "/work/mh0033/m300883/3rdPanel/data/indexDiffPattern/all_all_ind.nc"
    ).pc
    all_first_ind = xr.open_dataset(
        "/work/mh0033/m300883/3rdPanel/data/indexDiffPattern/all_first_ind.nc"
    ).pc
    all_last_ind = xr.open_dataset(
        "/work/mh0033/m300883/3rdPanel/data/indexDiffPattern/all_last_ind.nc"
    ).pc

    ind_fA, ind_fF, ind_fL, ind_lA, ind_lF, ind_lL = period_extreme(
        [all_all_ind, all_first_ind, all_last_ind]
    )


if __name__ == "__main__":
    main()
