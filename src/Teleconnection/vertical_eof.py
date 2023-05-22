import numpy as np
import xarray as xr
from tqdm.notebook import tqdm

import src.Teleconnection.rolling_eof as rolling_eof


def vertical_eof(
    xarr: xr.DataArray,  # the data to be decomposed
    nmode: int = 2,  # the number of mode to be generated
    window: int = 10,  # the window size if fix_pattern = 'False' is adopted.
    fixed_pattern: str = "all",  # the fixed_pattern mode
    independent: bool = True,  # the vertical strategy.
):
    """
    different way to do the eof vertically,
    **Arguments**:
        *xarr*: xr.DataArray,             # the data to be decomposed
        *nmode*: int = 2,                 # the number of mode to be generated
        *window*: int = 10,               # the window size if fix_pattern = 'False' is adopted.
        *fixed_pattern*: str= "all",      # the fixed_pattern mode
        *independent*: bool=True,         # the vertical strategy.
        *standard*: bool= True,           # standard pc or not
    **Return**

        *eof*, *pc* and *fra*
    """
    if independent == True:
        eof_result = independent_eof(
            xarr,
            nmode=nmode,
            window=window,
            fixed_pattern=fixed_pattern,
        )
    else:
        eof_result = dependent_eof(
            xarr,
            nmode=nmode,
            window=window,
            fixed_pattern=fixed_pattern,
        )

    return eof_result


def independent_eof(xarr, **kwargs):
    """
    do eof independently over all layers.
    **Arguments**
        *xarr* : the xarr to be composed.
        *fixed_pattern* : the method to generate pc.
        *method*: "rolling_eof" or "eof"
    **Return**
        EOF, PC and FRA.
    """

    print("     indenpendtly decomposing...")

    eof_result = xarr.groupby("plev").apply(
        rolling_eof.rolling_eof, **kwargs
    )

    return eof_result


def dependent_eof(xarr, **kwargs):
    """
    do eof independently over all layers.
    **Arguments**
        *xarr* : the xarr to be composed.
        *fixed_pattern* : the method to generate pc.
        *method*: "rolling_eof" or "eof"
    **Return**
        EOF, PC and FRA.
    """
    print("     dependently decomposign...")
    eof_result = rolling_eof.rolling_eof(xarr, **kwargs)

    return eof_result
