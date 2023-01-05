############ imports ###############
import numpy as np
import pandas as pd
import xarray as xr
from eofs.standard import Eof
from tqdm.notebook import tqdm, trange

import src.Teleconnection.tools as tools


def doeof(
    data: xr.DataArray,
    nmode: int = 2,
    dim: str = "com",
    standard: bool = True,
    shuffle: bool = True,
):
    """
    do eof to seasonal data along a combined dim

    **Arguments**:
        *seasondata*: The data to be decomposed, where the first dim should be the dim of 'com' or 'time'.
        *nmode*: how many modes, mode=2,means NAO and EA respectively.
        *dim*: along which dim to do the eof (calculate the co variance)
        *standard*: standard the output pc or not.
        *shuffle* random order the com (ens and time) or not.

    **Returns**:
        eof: DataArray, spatial patterns scaled (multiplied) with the temporal std of seasonal index.
             has the same spatial size as the input seasondata.[mode, lat,lon,...]
        pc: DataArray, seasonal index, if standard = True, scaled (divided) by the temporal std of itself. [mode,time]]
        exp_var: explained variance of each mode. [mode]

    """

    # make sure that the first dim is the 'com' or 'time'.
    try:
        data = data.transpose(dim, ...)
    except ValueError:
        print("no combined dimension found. use tools.stackens() first")

    # shuffle
    if shuffle:
        data = tools.random_order(data, dim=dim)

    # weights
    wgts = tools.sqrtcoslat(data)

    # EOF decompose
    solver = Eof(
        data.values, weights=wgts, center=True
    )  # if it's com dim, is is right to remove the mean
    # along the com dim?
    eof = solver.eofs(neofs=nmode)  # (mode,lat,lon,...)
    pc = solver.pcs(npcs=nmode)  # (com,mode)
    fra = solver.varianceFraction(nmode)  # (mode)

    # deweight
    eof_dw = eof / wgts

    # standarize coef
    std_pc = (np.std(pc, axis=0)).astype(
        "float64"
    )  # (mode)  # here should it be temporal std????
    dim_add_sp = np.hstack([nmode, tools.detect_spdim(data)])  # [-1,1,1] or [-1,1,1,1]
    std_pc_sp = std_pc.reshape(dim_add_sp)

    # eof should always be normalized.
    eof = eof_dw * std_pc_sp

    # while pc sometime should not for comparation.
    if standard:
        pc = pc / std_pc

    # xarray container for eof
    eof_cnt = data[:nmode]
    eof_cnt = eof_cnt.rename({dim: "mode"})
    eof_cnt["mode"] = ["NAO", "EA"]

    # to xarray
    eofx = eof_cnt.copy(data=eof)
    pcx = xr.DataArray(
        pc, dims=[dim, "mode"], coords={dim: data[dim], "mode": ["NAO", "EA"]}
    )
    frax = xr.DataArray(fra, dims=["mode"], coords={"mode": ["NAO", "EA"]})
    eofx.name = "eof"
    pcx.name = "pc"
    frax.name = "exp_var"

    # change sign
    coef = sign_coef(eofx)
    eofx = eofx * coef
    pcx = pcx * coef

    # unstack the dim 'ens' and 'time' or 'win'
    pcx = pcx.unstack()

    return eofx, pcx, frax


def project_field(fieldx, eofx, dim="com", standard=True):
    """project original field onto eofs to get the temporal index.

    Different from python eofs package, here if there are three dimensions in sptial,
    i.e, [lat,lon,height], the projected pc is calculated independently from each height.

    **Arguments:**

        *field*: the DataArray field to be projected
        *eof*: the eofs
        *standard*: whether standardize the ppc with its std or not

    **Returns:**

        projected pcs
    """
    fieldx = fieldx.transpose(dim, ...)
    neofs = eofx.shape[0]

    # weight
    wgts = tools.sqrtcoslat(fieldx)
    field = fieldx.values * wgts

    # fill with nan
    try:
        field = field.filled(fill_value=np.nan)
    except AttributeError:
        pass

    # flat field to [time,lon-lat] or [time,lon-lat,heith]
    records = field.shape[0]
    channels = np.product(field.shape[1:3])  # only lat and lon stack here.
    nspdim = len(tools.detect_spdim(fieldx))  # how many spatial dims
    if nspdim > 2:
        heights = eofx.shape[3]
    try:
        field_flat = field.reshape([records, channels, heights])
    except NameError:
        field_flat = field.reshape([records, channels])

    # non missing value check
    nonMissingIndex = np.where(np.logical_not(np.isnan(field_flat[0])))[0]
    field_flat = field_flat[:, nonMissingIndex]

    # flat eof to [mode, space]
    try:
        _flatE = eofx.values.reshape(neofs, channels, heights)
    except NameError:
        _flatE = eofx.values.reshape(neofs, channels)

    eofNonMissingIndex = np.where(np.logical_not(np.isnan(_flatE[0])))[0]

    # missing value align check
    if (
        eofNonMissingIndex.shape != nonMissingIndex.shape
        or (eofNonMissingIndex != nonMissingIndex).any()
    ):
        raise ValueError("field and EOFs have different " "missing value locations")
    eofs_flat = _flatE[:, eofNonMissingIndex]

    # for three dimentional space data
    try:
        projected_pcs = []  # for all height layers
        for h in range(heights):
            field_flat_h = field_flat[:, :, h]
            eofs_flat_h = eofs_flat[:, :, h]
            projected_pc = np.dot(field_flat_h, eofs_flat_h.T)
            projected_pcs.append(projected_pc)
        projected_pcs = np.array(projected_pcs)

        PPC = xr.DataArray(
            projected_pcs,
            dims=[
                fieldx.dims[-1],
                fieldx.dims[0],
                eofx.dims[0],
            ],  # [height,record,mode]
            coords={
                fieldx.dims[-1]: fieldx[fieldx.dims[-1]],
                fieldx.dims[0]: fieldx[fieldx.dims[0]],
                eofx.dims[0]: eofx[eofx.dims[0]],
            },
        )
    # for 2-d space
    except NameError:
        projected_pcs = np.dot(field_flat, eofs_flat.T)
        PPC = xr.DataArray(
            projected_pcs,
            dims=[fieldx.dims[0], eofx.dims[0]],
            coords={
                fieldx.dims[0]: fieldx[fieldx.dims[0]],
                eofx.dims[0]: eofx[eofx.dims[0]],
            },
        )
    PPC.name = "pc"

    # to unstack 'com' to 'time' and 'ens' if 'com' exists.
    PPC = PPC.unstack()

    if standard:
        PPC = tools.standardize(PPC)
    return PPC


def sign_coef(eof):
    """
    function to calculate the coefficient for eof, so that the sign is consistent.
    for NAO, the positive NAO with a low in the North and high at the south.
    for EA, the positive EA with a low in the center.
    **Arguments**:
        eof: the eof (spatial pattern) to be changed. much include NAO and EA the same time.
    **Returns**:
        coefficient of NAO and EA in xarray.
    """
    # NAO
    coef_NAO = (
        eof.sel(lat=slice(90, 60), lon=slice(-70, -10), mode="NAO").mean(
            dim=["lat", "lon"]
        )
        < 0
    )
    coef_NAO = 2 * coef_NAO - 1  # to make 1 to 1 , 0 to -1

    # EA
    coef_EA = (
        eof.sel(lat=slice(65, 45), lon=slice(-40, 40), mode="EA").mean(
            dim=["lat", "lon"]
        )
        < 0
    )
    coef_EA = 2 * coef_EA - 1

    return xr.concat([coef_NAO, coef_EA], dim="mode")
