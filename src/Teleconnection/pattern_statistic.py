import numpy as np
import pandas as pd


def select_co2(xarr):
    """
    select 3 CO_2 stages:
        1856: $1 \times CO_2$
        1921: $2 \times CO_2$
        1991: $4 \times CO_2$
    **Arguments**
        *xarr* the DataArray to be selected.
    **Return**
        *xarr* the DataArray only with three years.
    """
    co2time = ["1856-03-03", "1921-03-16", "1991-03-16"]
    co2time = pd.DatetimeIndex(co2time)
    xarr_CO2 = xarr.sel(time=xarr.time.dt.year.isin(co2time.year))
    return xarr_CO2


def lon_height(eof, mode="EA"):
    """
    calculate zonally mean of eof
    **Arguments**
        *eof* The eofs to be averaged [time,lat,lon,plev]
    **Return**
        eof as a function of longitude and height
    """
    # height
    eof_height = eof.sel(hlayers=slice(20000, 100000))
    eof_height["hlayers"] = eof_height["hlayers"] / 100  # from pa to hpa

    # lats
    eof_lat = eof_height.sel(lat=slice(80, 30))

    # groupby
    lon_bins = np.arange(-90, 41, 5)
    EA_lon_height = (
        eof_lat.sel(mode=mode).groupby_bins("lon", bins=lon_bins).mean(dim="lat")
    )
    return EA_lon_height


def lat_height(eof, mode="NAO"):
    """
    calculate meridional mean of eof
    **Arguments**
        *eof* the eofs to be averaged.
        *mode* which mode to be calculated
    **Return**
        eof as function of lattitude and height.
    """

    # height
    eof_height = eof.sel(hlayers=slice(20000, 100000))
    eof_height["hlayers"] = eof_height["hlayers"] / 100

    # groupby
    lat_bins = np.arange(20, 81, 4)
    lat_labels = np.arange(22, 81, 4)
    EA_lat_height = (
        eof_height.sel(mode=mode)
        .groupby_bins("lat", bins=lat_bins, labels=lat_labels)
        .mean(dim="lon")
    )
    return EA_lat_height
