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


def mode_return_period(index: xr.DataArray, mode: str, hlayers: int = 50000):
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

    # pos
    all_pos = return_period(index, "pos")
    first10_all_pos, last10_all_pos = tenYear_return_period(all_pos)
    first10_median_pos, last10_median_pos = tenYear_return_period_median(
        all_pos, method="max"
    )

    # neg
    all_neg = return_period(index, "neg")
    first10_all_neg, last10_all_neg = tenYear_return_period(all_neg)
    first10_median_neg, last10_median_neg = tenYear_return_period_median(
        all_neg, method="min"
    )

    return (
        first10_all_pos,
        last10_all_pos,
        first10_median_pos,
        last10_median_pos,
        first10_all_neg,
        last10_all_neg,
        first10_median_neg,
        last10_median_neg,
    )


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
        (
            _,
            _,
            POS_median_first,
            POS_median_last,
            _,
            _,
            NEG_median_first,
            NEG_median_last,
        ) = mode_return_period(index, mode, hlayers)

        pos[i, 0] = POS_median_first["return period"].values
        pos[i, 1] = POS_median_last["return period"].values

        neg[i, 0] = NEG_median_first["return period"].values
        neg[i, 1] = NEG_median_last["return period"].values
    return pos, neg


def median_return_period(d, method):
    """
    the median return period
    """
    ranks = d["return period"].rank(pct=True, method=method)
    close_to_median = abs(ranks - 0.5)
    return d.loc[[close_to_median.idxmin()], :]


def tenYear_return_period(return_periods):
    """
    split into first10 and last10
    """
    first10 = return_periods.loc["1850":"1860"]
    last10 = return_periods.loc["1990":"1999"]
    return first10, last10


def tenYear_return_period_median(all_return_period, method):
    """
    calculate the median of first10 and last10.
    """
    first10, last10 = tenYear_return_period(all_return_period)
    first10_median = median_return_period(first10, method=method)
    last10_median = median_return_period(last10, method=method)
    return first10_median, last10_median
