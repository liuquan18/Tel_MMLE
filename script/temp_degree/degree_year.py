import xarray as xr
import numpy as np

def return_year(xarr):
    return xarr.time.dt.year.values

def degree_year(fldmean:xr.DataArray):
    """
    to calculate the year when the mean global surface temperature reaches 1,2,and 4 degrees.
    **Argument**
        *fldmean* the fldmean of tsurf
    """
    if isinstance(fldmean,xr.DataArray):
        pass
    else:
        print("only DataArray is accapted, DataSet recevied")

    try:
        fldmean.lon.size == 1 & fldmean.lat.size == 1
    except ValueError:
        print("the fldmean temperature should be calculated first")

    # ens mean
    if fldmean.ens.size != 1:
        fld_ens_mean = fldmean.mean(dim = 'ens')
    else:
        fld_ens_mean = fldmean

    # squeeze
    mean = fld_ens_mean.squeeze()

    # anomaly
    anomaly = mean-mean[0]

    years = np.zeros(3)

    # 0 degree (1855)
    years[0] = return_year(anomaly[5])

    # 2 degree
    years[1] = return_year(anomaly.where(anomaly>=2,drop=True)).squeeze()[0]

    # 4 degree
    years[2] = return_year(anomaly.where(anomaly>=4,drop=True)).squeeze()[0]

    return years