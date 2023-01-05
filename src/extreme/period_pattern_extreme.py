import xarray as xr
import src.extreme.extreme as ext

import importlib

importlib.reload(ext)


def normalize(period_period: xr.DataArray, all_all: xr.DataArray,dim = 'time'):
    """
    normalize the temporal index of first10 and last10 years
    with the mean and std of all index from all pattern
    to make the index in the first10 years comparable with the
    index in the last10 years.
    **Arguments**
        *period_period* the index of the first10 or last10 years
        *all_all* the index of all the years with all pattern
        *dim* along which dim to do the normalize.
    **Return**
        *peiod_period* the normalized index
    """
    mean = all_all.mean(dim= dim)
    std = all_all.std(dim=dim)
    period_period = (period_period - mean) / std
    return period_period


def period_extreme(peiod_index: xr.DataArray, threshold=2):
    """
    mask out non-extreme cases as np.nan.
    **Arguments**
        *period_index* the three index of all years from the three patterns.
                    should be order in [all_all, all_first, all_last]
    **Return**
        *period_extreme* with all-non-extreme labeled as np.nan
    """
    period_extreme = ext.extreme(peiod_index, threshold)
    return period_extreme


def _period_extreme_count(period_extreme: xr.DataArray, dim: tuple):
    """
    count the number of the extremes
    **Arguments**
        *period_extrme* the index with non-extreme masked as np.nan from period_extreme function
        *dim* along which dim to count the extremes
    **Return**
        *extreme_count* the count of the extremes
    """
    return ext.count_extreme(period_extreme, dim=dim)


def period_extreme_count(
    period_period: xr.DataArray,
    all_all: xr.DataArray = None,
    dim=("time", "ens"),
    threshold=2,
    standard: str = True,
):
    """
    the extreme count of ten_year period
    **Arguments**
        *period_period* the index of the first10-first or last10-last
        *all_all* the index of all-all pattern
        *dim* along which dim to count.
        *threshod* threshold for extreme
        *standard*
    **Return**
        *extreme_count* the count of the extremes in this ten years.
    """
    if all_all is not None:
        if standard:
            period_period = normalize(period_period, all_all)

    extreme = period_extreme(period_period, threshold)
    count = _period_extreme_count(extreme, dim=dim)
    return count
