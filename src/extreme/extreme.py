"""
extreme pacakge
"""
import xarray as xr


def extreme(xarr, threshod=2):
    """
    mask out non-extreme data.
    label non-extreme data as np.nan.
    **Arguments**
        *xarr* the xarr to be process
    **Return**
        *ex* extreme xarray with one extra dimension called 'extr_type',non
        extreme data labeled as np.nan
    """
    pos_ex = xarr.where(xarr > threshod)
    neg_ex = xarr.where(xarr < -1 * threshod)
    ex = xr.concat([pos_ex, neg_ex], dim=["pos", "neg"])
    ex = ex.rename({"concat_dim": "extr_type"})
    return ex


def count_extreme(extreme_nan, dim=("time", "ens")):
    """
    count the number of extreme cases in xarr
    **Arguments**
        *extreme_nan* the xarr where the non-extreme points are labeled as np.nan
    **Return**
        number of extreme cases, with the coordinate of 'hlayer' and 'mode'
        reserved.
    """
    count = extreme_nan.count(dim=dim)

    return count

