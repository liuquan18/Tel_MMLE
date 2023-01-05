import numpy as np
import xarray as xr
from tqdm.notebook import tqdm

import src.Teleconnection.rolling_eof as rolling_eof


def vertical_eof(
        xarr: xr.DataArray,             # the data to be decomposed
        nmode: int = 2,                 # the number of mode to be generated
        window: int = 10,               # the window size if fix_pattern = 'False' is adopted.
        fixed_pattern: str= "all",      # the fixed_pattern mode
        independent: bool=True,         # the vertical strategy. 
        standard: bool= True,           # standard pc or not
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
        eof, pc, fra = independent_eof(
            xarr, 
            nmode = nmode, 
            window = window, 
            fixed_pattern = fixed_pattern,
            standard = standard)
    else:
        eof, pc, fra = dependent_eof(
            xarr, 
            nmode = nmode, 
            window = window, 
            fixed_pattern = fixed_pattern,
            standard = standard)

    return eof, pc, fra

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
    eofs = []
    pcs = []
    fras = []

    hlayers = xarr.hlayers
    for h in tqdm(hlayers):
        field = xarr.sel(hlayers=h)
        eof, pc, fra = rolling_eof.rolling_eof(field, **kwargs)

        eofs.append(eof)
        pcs.append(pc)
        fras.append(fra)
    eofs = xr.concat(eofs, dim=xarr.hlayers)
    pcs = xr.concat(pcs, dim=xarr.hlayers)
    fras = xr.concat(fras, dim=xarr.hlayers)
    return eofs, pcs, fras


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
    eofs, pcs, fras = rolling_eof.rolling_eof(xarr, **kwargs)

    return eofs, pcs, fras

