import numpy as np
import pandas as pd
import xarray as xr
from tqdm.notebook import tqdm, trange

import src.Teleconnection.spatial_pattern as ssp
import src.Teleconnection.tools as tools


def rolling_eof(xarr, nmode=2, window=10, fixed_pattern="all", standard=True):
    """do eof analysis with in a rolling window.

    rolling EOF is like rolling mean, here the default window is 10 years. The EOFs, PCs and
    exp_var of one specifice year is decomposed from a combined dataset,which is generated by
    combing the data 5 years before and 4 years later (totaly 10 years of data), and by concating
    the ensemble dim together, we have a series of 1000 length.

    **Arguments**:

        *xarr*: the DataArray that is going to be decomposed by rolling_eof. [time,lat,lon,ens,height,decade]
        *nmode*: the number of modes reserved.
        *window*: the rolling window.
        *fixed_pattern*: string, determine how the pcs generated.
                       - if fixed_pattern = 'all', then 'pc' is calculated by projecting the original
                       fields onto a temporally-fixed pattern. i.e, the 'pc' is the second output
                       (pcx) of 'doeof' on a combined dataset, which concates all the  ensembles
                       and years togerther. In which case the spatial pattern is assumed to be
                       fixed.
                       - if fixed_pattern = 'first'. then 'pc' is calculted by projecting the original
                       fields onto the spatial pattern of the first 10 years.
                       - if fixed_pattern = 'last', the 'pc' is calculated by projecting the original
                       fields onto the spatial pattern of the last 10 years.
                       - if fixed_pattern = False, the pc is calculated by projecting the seasonal
                       data onto to the spatail pattern of this year, in this case, the sptial
                       pattern is gotten by the data ten years around this year. And the length
                       of the index should be shorter than the original series, since not enough
                       years in the begining and the end.
        *return_full_eof*: whether to return the full rolling eofs. if False, only the first and the
                           last eof is returned. if True, it would take a long time.
        *standard* whether the pc (from project field) should be standard with its own std and mean
    **Returns**:

        EOF: The eofs, but now with the first dim as the time.[time-10,mode,lat,lon,(height)]
        PC: the pcs, if fixed_pattern=True, the pcs should be [ens,time,mode]
                     if fixed_pattern=False, the pcs should be [ens,time-10, mode]
        FRA: the explained variances, the shape should be [time-10,mode]
    """

    # the validtime period where totally ten years of data are fully avaiable.
    gap = int(window / 2)
    validtime = xarr.isel(time=slice(gap, -1 * gap)).time

    # if do the all-all decompose
    if fixed_pattern == "all":  # a little different from the following two.
        print("     using the all pattern")
        eof_result = ssp.doeof(
            tools.stack_ens(xarr, withdim="time"),
            nmode=nmode,
            dim="com",
            standard=False,
        )
        # here the pc is not directly used since the eof is multiplied by the std of pc, then if we
        # do the project-field, the resulted projectd-pc is not the same as the pc from the solver.
        # in order to make it the same  order as the following, we do project-field to get the index.
        PC = fixed_pc(xarr, eof_result["eof"], standard=standard)

    elif fixed_pattern == "first":
        # only the EOF of the first10 is needed.
        print("     using the first pattern")
        eof_result = changing_eofs(xarr, validtime[0], nmode=nmode, window=window)
        PC = fixed_pc(xarr, eof_result["eof"], standard=standard)

    elif fixed_pattern == "last":
        print("     using the last pattern")
        eof_result = changing_eofs(xarr, validtime[-1], nmode=nmode, window=window)
        PC = fixed_pc(xarr, eof_result["eof"], standard=standard)

    elif fixed_pattern == "decade":
        print("     decomposing everty ten years")
        # select the middle year of each decade from validtime
        decade_time = validtime[::window]
        eof_result = changing_eofs(xarr, decade_time, nmode=nmode, window=window)
        PC = decadal_pc(xarr, decade_time, eof_result["eof"],window=window)

    elif fixed_pattern == "False":
        print("     no fixed pattern used")
        eof_result = changing_eofs(xarr, validtime, nmode=nmode, window=window)
        PC = changing_pc(xarr, validtime, eof_result["eof"], standard=standard)

    # replace the pc with the pc from projection
    eof_result["pc"] = PC
    return eof_result


def changing_eofs(xarr, validtime, nmode, window):
    """getting the rolling eofs and fras from xarr.

    **Arguments**

        *xarr* : the DataArray to be composed.
        *validtime*: the times who has the full 10 years periods around it.
        *nmode* : how many modes to decompose.
        *windo* : rolling window.

    **Return**

        eof and fra.
    """

    field = rolling(xarr, win=window)
    field = field.stack(com=("ens", "decade"))
    # select only the valid time
    field = field.sel(time=validtime)

    # changing eofs (dynamic):
    if validtime.size > 1:
        eof_result = field.groupby("time").apply(
            lambda x: ssp.doeof(x, nmode=nmode, dim="com")
        )

    # for one pattern (first,all,last) and all time step
    elif validtime.size == 1:
        eof_result = ssp.doeof(field, nmode=nmode, dim="com")

    # drop vars
    try:
        eof_result = eof_result.drop_vars("decade")
    except:
        pass

    return eof_result


def fixed_pc(xarr, pattern, standard, dim="time"):
    """projecting the xarr to a fixed spatial pattern.

    **Arguments**

        *xarr*: the xarr to be decomposed.
        *pattern*: the fixed pattern to project on.

    **Returns**

        *pcx*: the coresponding temporal index
    """
    # stack
    fieldx = tools.stack_ens(xarr, withdim=dim)
    pc = ssp.project_field(fieldx, pattern, standard=standard)
    return pc


def decadal_pc(xarr, validtime, EOF, window):
    """decompose the xarr into decadal patterns.

    **Arguments**

        *xarr* : the DataArray to be composed.
        *nmode* : how many modes to decompose.
        *windo* : rolling window.

    **Return**

        eof and fra.
    """

    field = rolling(xarr, win=window)
    PC = []
    # select only the valid time
    for time in validtime:
        lowyear = pd.Timestamp(time.values) - pd.DateOffset(years = window/2)
        highyear = pd.Timestamp(time.values) + pd.DateOffset(years = window/2 - 1)
        time_window = pd.date_range(lowyear, periods=window,freq='A-MAR')

        field_dec = field.sel(time=time)
        pattern = EOF.sel(time=time)
        pc = fixed_pc(field_dec, pattern, standard=False, dim="decade")

        # change 'decade' to 'time'
        pc = pc.drop_vars("time")
        pc = pc.rename({"decade": "time"})
        pc['time'] = time_window
        PC.append(pc)
    PC = xr.concat(PC, dim='time')
    return PC


def changing_pc(xarr, validtime, EOF, standard):
    """projecting the xarr to a changing spatial pattern.
    The patterns can be a
        pattern for each year (dim = 'ens')
        pattern for each decade (dim = 'decade')

    **Arguments**

        *xarr* the array to be composed.
        *validatime* the timeseries who has full 10 years around it.
        *EOF* the changing EOFs.

    **Return**

        changing pc, whose length is time-10.
    """
    PC = []
    for time in validtime:
        field = xarr.sel(time=time)
        pattern = EOF.sel(time=time)

        # here the standard must be False since there is only one time step here.
        pc = ssp.project_field(
            field,
            pattern,
            dim="ens",
            standard=False,
        )  # project all ens onto one eof.
        PC.append(pc)
    PC = xr.concat(PC, dim=validtime)

    # The real standardization is done here.
    if standard:
        PC = tools.standardize(PC)

    return PC


def rolling(xarr, win=10):
    """rolling the xarr with a time window of 'win'.

    **Arguments**:

        xarr: DataArray to be rolled.
        win: the window size

    **Return**:

        DataArray being rolled.
    """
    # using the rolling func in xarray.
    roller = xarr.rolling(time=win, center=True)

    # construct the rolling object to a DatArray
    rolled = roller.construct("decade")

    return rolled
