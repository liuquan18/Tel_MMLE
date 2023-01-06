from pyextremes import get_extremes, get_return_periods
import xarray as xr
import numpy as np


def return_period(index: xr.DataArray, extreme_type: str):
    """
    return return period given the extreme_type with EVT
    **Parameters**:
        *index* the index for all years and all ens
        *extreme_type* the extreme type: pos or nega
    **Return**:
        *returnPeriod* the return period from BM method.
    """
    # get the extremes

    if extreme_type == "pos":
        extremes_type = "high"
        extremes = index.max(dim="ens").to_dataframe()["pc"]

    elif extreme_type == "neg":
        extremes_type = "low"
        extremes = index.min(dim="ens").to_dataframe()["pc"]

    returnPeriod = get_return_periods(
        ts=extremes,
        extremes=extremes,
        block_size="365.2425D",
        return_period_size="365.2425D",
        plotting_position="cunnane",
        extremes_type=extremes_type,
        extremes_method="BM",
    )
    return returnPeriod


def mode_return_period(
    index: xr.DataArray, mode: str, periods: list, hlayers: int = 50000
):
    """
    the return period of {mode} at one altitude layer.
    split to first10 and last10 years. calculate the media return period.
    **Parameters**
        *index* index of NAO and EA in xr.DataArray.
        *mode* the mode NAO or EA.
        *hlayers* the height layer to be calculate
    **Return**
        pos and neg for first10, last10, and media return period. (8 output...)
    """

    index = index.sel(mode=mode, hlayers=hlayers)

    period_pos, media_period_pos = split_period(index, periods,'pos')
    period_neg, media_period_neg = split_period(index,periods, 'neg')

    return (period_pos, media_period_pos, period_neg, media_period_neg)

def split_period(index, periods,extr_type):
    if extr_type == 'pos':
        method = 'max'
    elif extr_type == 'neg':
        method = 'min'
    all = return_period(index, extr_type)
    period = [all.loc[period] for period in periods]
    media_period = [median_return_period(p, method) for p in period]
    return period,media_period


def median_return_period(d, method):
    """
    the median return period
    """
    ranks = d["return period"].rank(pct=True, method=method)
    close_to_median = abs(ranks - 0.5)
    return d.loc[[close_to_median.idxmin()], :]


def vertical_return_period(index: xr.DataArray, mode: str):
    """
    the media return period of all altitudes.
    **Parameters**:
        *index* the index of NAO and EA.
        *mode* NAO or EA.
    **Return**
        *pos* the vertical profile of media return period for positive extremes.
        *neg* the vertical profile of medai return period for negative extremes.
    """

    pos = np.zeros((index.hlayers.size, 2))
    neg = np.zeros((index.hlayers.size, 2))

    for i, hlayers in enumerate(index.hlayers):
        _,pos[i],_,neg[i]= mode_return_period(index, mode, hlayers)
    return pos, neg

