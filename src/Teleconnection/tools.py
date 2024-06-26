import numpy as np
import pandas as pd
import xarray as xr
from sklearn.utils import shuffle


########### functions to do pre-process ###########
def split_ens(xarr):
    """split a dataset generated by cdo -apply, where the 'plev' and 'ens' are mixed into 'plev' dim.

    **Arguments**:

        xarr: DataArray with plev and ens mixed.

    **Returns**:

        DataArray with splited dims.

    """

    # one dim: plev, i.e, the vetical gph of MPI-GE, the ture plev coordinate.
    hlayers = np.array(
        [
            1.00e05,
            9.25e04,
            8.50e04,
            7.75e04,
            7.00e04,
            6.00e04,
            5.00e04,
            4.00e04,
            3.00e04,
            2.50e04,
            2.00e04,
            1.50e04,
            1.00e04,
            7.00e03,
            5.00e03,
            3.00e03,
            2.00e03,
            1.00e03,
            7.00e02,
            5.00e02,
            3.00e02,
            2.00e02,
            1.00e02,
            5.00e01,
            2.00e01,
            1.00e01,
        ]
    )

    # The other dim: ens, i.e, the index for ens
    ens = np.arange(100)

    # createa a multiindex which has the same length as the mix dim in xarr.
    ind = pd.MultiIndex.from_product((ens, hlayers), names=("ens", "hlayers"))

    # replace the old coords with ind, and then unstack.
    replace = xarr.assign_coords(plev=ind).unstack("plev")
    replace = replace.rename({'hlayers':'plev'})

    return replace


def stack_ens(xarr, withdim="decade"):
    """
    The first dim of input data for python-eofs-package standard interface should be 'time', but
    here we do eof not along the time dim only, but the win (10) years of all ensembles. so should
    stack the dims 'win' and 'ens' together first.
    **Arguments**:
        xarr: the rolled DataArray to be stacked.
        withdim: with which dim the ens should be stacked with.
    Return:
        xarr: with the dim 'win' or 'time' and 'ens' combined to 'com'.
    """

    time_com_space = xarr.stack(com=("ens", withdim))
    return time_com_space


def standardize(xarr,dim = None):
    """
    standardize the DataArray with the temporal mean and std. Note here the function standardize
    the input with the mean and std of its own. for comparation the mean and std should be the
    same across different enterties.
    or the function here should be applied onto the input, not the output.
    **Arguments**:
        xarr: The DataArray to be standarized.
    **Returns**:
        xarr: standarized DataArray
    """
    if dim is None:
        try:
            time_mean = xarr.mean(dim="time")
            time_std = xarr.std(dim="time")
        except ValueError:
            time_mean = xarr.mean(dim="decade")
            time_std = xarr.std(dim="decade")
    else:
        time_mean = xarr.mean(dim = dim)
        time_std = xarr.std(dim = dim)
    
    return (xarr - time_mean) / time_std


########## Function to do EOF #####################

def sqrtcoslat(xarr):
    """
    calculte the square-root of the cosine of the latitude as the weight
    **Arguments**:
        xarr: the DataArray to calculate the weight, with the shape ['com','lat','lon',...]
    **Returns**:
        weight with the right shape.
    """
    # weight values
    wgts = np.sqrt(np.cos(np.deg2rad(xarr.lat)))
    
    # make the shape of wgts the same as xarr
    W = xr.ones_like(xarr)
    wgts = W * wgts
    return wgts


 