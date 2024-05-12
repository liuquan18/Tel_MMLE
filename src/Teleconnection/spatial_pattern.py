############ imports ###############
import numpy as np
import xarray as xr
from eofs.standard import Eof

import src.Teleconnection.tools as tools


def doeof(
    data: xr.DataArray,
    nmode: int = 2,
    dim: str = "com",
    standard: str = 'eof_spatial_std',
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

    # make sure that there are only three dims
    data = data.squeeze()

    # make sure that the first dim is the 'com' or 'time'.
    try:
        data = data.transpose(dim, ...)
    except ValueError:
        print("no combined dimension found. use tools.stackens() first")

    # weights
    wgts = tools.sqrtcoslat(data) # same coords as data

    # EOF decompose
    solver = Eof(
        data.values, weights=wgts, center=False
    )  
    
    eof = solver.eofs(neofs=nmode)  # (mode,lat,lon,...)
    pc = solver.pcs(npcs=nmode)  # (com,mode)
    fra = solver.varianceFraction(nmode)  # (mode)

    # eof to xarray
    eofx, pcx, frax = eofs_to_xarray(data, eof, pc, fra)

    # deweight
    eofx = eofx / wgts[0]

    # standardize, here the loading gives to the pc, to make the index from different spatil pattern comparable.
    if standard == "eof_spatial_std":
        eofx, pcx = standard_by_eof_spatial_std(eofx, pcx)
    elif standard == "pc_temporal_std":
        eofx, pcx = standard_by_pc_temporal_std(eofx, pcx)

    # fix the sign, so that the North center of action is always low.
    if eofx.mode.size == 2:
        eofx, pcx = fix_sign(eofx, pcx) # only when the first two modes are decomposed

    # make sure at the loc where the data==np.nan, the eof==np.nan as well.
    map_data = data[0]  # just one map
    eofx = eofx.where(np.logical_not(map_data.isnull()), map_data)

    # unstack the dim 'ens' and 'time' or 'win'
    pcx = pcx.unstack()

    # dorp vars
    try:
        eofx = eofx.drop_vars(("ens", "time", dim))
    except ValueError:
        pass

    # to dataset
    eof_result = xr.Dataset({"eof": eofx, "pc": pcx, "fra": frax})

    return eof_result

def standard_by_eof_spatial_std(eofx, pcx):
    std_eof = eofx.std(dim=("lat", "lon"))
    eofx = eofx / std_eof
    std_eof = std_eof.squeeze()
    pcx = pcx * std_eof
    return eofx,pcx

def standard_by_pc_temporal_std(eofx, pcx):
    std_pc = pcx.std(axis = 0) # either along 'com' or 'time'
    pcx = pcx / std_pc
    std_pc = std_pc.squeeze()
    eofx = eofx * std_pc
    return eofx,pcx

def eofs_to_xarray(data, eof, pc, fra):
    if pc.shape[1] ==2:
        modes = ['NAO','EA']
    else:
        modes = np.arange(pc.shape[1])

    reduce_dim = data.dims[0] # 'com' or 'time'
    eof_cnt = data[0]
    time_tag = data.unstack().time.values[0] # the time tag of the first map
    eof_cnt = eof_cnt.drop_vars(reduce_dim) # drop the dim 'com' or 'time'
    eof_cnt = eof_cnt.expand_dims(dim = {'mode':modes,'decade':[time_tag]},axis = [0,-1]) # add the dim 'mode' and 'decade'

    eof = eof[..., np.newaxis] # add one new dimension for the info of decade (time of the beginning of the decade)
    fra = fra[..., np.newaxis] # add one new dimension for the info of decade (time of the beginning of the decade)
    eofx = eof_cnt.copy(data=eof)

    # pc to xarray
    pcx = xr.DataArray(
        pc, dims=[reduce_dim, "mode"], coords={reduce_dim: data[reduce_dim], "mode":modes}
    )
    frax = xr.DataArray(fra, dims=["mode","decade"], coords={"mode": modes,'decade':[time_tag]})
    return eofx,pcx,frax


def fix_sign(eof,pc):
    """
    function to calculate the coefficient for eof, so that the sign is consistent.
    for NAO, the positive NAO with a low in the North and high at the south.
    for EA, the positive EA with a low in the center.
    **Arguments**:
        eof: the eof (spatial pattern) to be changed. much include NAO and EA the same time.
    **Returns**:
        coefficient of NAO and EA in xarray.
    """

    # sortby lat since some dataset the lat goes from higher to lower.

    # NAO
    # check if the lat of the eof is from lower to higher
    NAO_box = neg_center_box(eof.sel(mode = 'NAO'), 60, 75, -75, -50)
    coef_NAO = (NAO_box.mean(dim=["lat", "lon"])< 0)
    coef_NAO = 2 * coef_NAO - 1  # to make 1 to 1 , 0 to -1

    # EA
    EA_box = neg_center_box(eof.sel(mode = 'EA'), 45, 65, -40, 40)
    coef_EA = (EA_box.mean(dim=["lat", "lon"])< 0)
    coef_EA = 2 * coef_EA - 1

    # change sign
    coef = xr.concat([coef_NAO, coef_EA], dim="mode")
    eof = eof * coef
    coef = coef.squeeze()
    pc = pc * coef
    return eof, pc

def neg_center_box(xarr, blat, tlat, llon, rlon): 
    """
    check if the lat and lon of the eof is from lower to higher.
    input the order of expected box in # bottom lat, top lat, left lon, right lon for the box
    return the order of the box in the eof.
    """
    if xarr.lat[0] > xarr.lat[-1]: # descending
        start_lat = tlat
        end_lat = blat
    else: # ascending
        start_lat = blat
        end_lat = tlat

    if xarr.lon[0] > xarr.lon[-1]:
        start_lon = rlon
        end_lon = llon

    else:
        start_lon = llon
        end_lon = rlon

    # if the longitude from 0-360, change the lon to -180-180
    if xarr.lon.min() >= 0 and xarr.lon.max() <= 360:
        if start_lon < 0:
            start_lon = start_lon + 360
        if end_lon < 0:
            end_lon = end_lon + 360

    
    return xarr.sel(lat=slice(start_lat, end_lat), lon=slice(start_lon, end_lon))

def project_field(fieldx, eofx, dim="com", standard=False):
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
    eofx = eofx.squeeze() # remove the 'decade' dim if exists

    # weight
    wgts_f = tools.sqrtcoslat(fieldx)
    field = fieldx * wgts_f

    wgts_eof = tools.sqrtcoslat(eofx)
    eofx = eofx * wgts_eof

    # fill with nan
    try:
        field = field.filled(fill_value=np.nan)
    except AttributeError:
        pass

    # flat field to [time,lon-lat] or [time,lon-lat,heith]
    field_flat = field.stack(spatial = ('lon','lat'))

    eof_flat = eofx.stack(spatial = ('lon','lat'))

    projected_pcs = np.dot(field_flat, eof_flat.T)
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