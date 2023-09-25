##This module contains functions that calculate blocking indices from xarray dataarrays of geopotential##

#%%
import numpy as np
import xarray as xr
import pandas as pd
#%%

##1D Tibaldi Molteni Index############

#get geopotential gradient in South region
def _GHGS(Z,lat_n=80,lat_0=60,lat_s=40,delta=0,MGI=False):
    
    Z_0=Z.sel(lat=lat_0+delta,method='nearest')
    Z_s=Z.sel(lat=lat_s+delta,method='nearest')
    ghgs=(Z_0-Z_s)/(lat_0-lat_s)
    return ghgs

#get geopotential gradient in North region
def _GHGN(Z,lat_n=80,lat_0=60,lat_s=40,delta=0):
    
    Z_0=Z.sel(lat=lat_0+delta,method='nearest')
    Z_n=Z.sel(lat=lat_n+delta,method='nearest')
    ghgn=(Z_n-Z_0)/(lat_n-lat_0)
    return ghgn

#Calculate the 1D Tibaldi Molteni blocking index:
def TM_index(Z, lat_n=80, lat_0=60, lat_s=40, deltas=[-5, 0, 5], thresh=-10):
    
    # calculate the mean over latitude
    TM = Z.mean(dim='lat')
    
    # set up memory storage
    s_test = np.zeros([len(deltas), *TM.shape])
    n_test = np.zeros([len(deltas), *TM.shape])

    # for each delta check if conditions are true
    for i, d in enumerate(deltas):
        
        # 1.) is the Southern gradient positive?
        ghgs = _GHGS(Z, lat_n, lat_0, lat_s, d)
        s_test[i] = ghgs.values > 0
        
        # 2.) is the Northern gradient less than -10m/deg lat?
        ghgn = _GHGN(Z, lat_n, lat_0, lat_s, d)
        n_test[i] = ghgn.values < thresh
        
    # if both conditions are true than the test is passed:
    test = np.array(s_test) * np.array(n_test)
    
    # if the test is passed for any delta we are happy:
    test = np.any(test, axis=0)
    
    # pack the result in an xarray DataArray and return
    TM = xr.DataArray(test, coords=Z.coords, dims=Z.dims)
    TM.name = "TM index"
    TM.attrs['long_name'] = "TM index"
    return TM

## 2D indices like Davini 2012##################

#Used to filter out low latitude blocks (LLB)
def _GHGS2(Z, lat_n=80, lat_0=60, lat_s=40, delta=0):
    
    Z_ss = Z.sel(lat=lat_s-15,method = 'nearest')
    Z_s = Z.sel(lat=lat_s,method = 'nearest')
    ghgs2 = (Z_s - Z_ss) / 15
    return ghgs2

#Calculate the 2D instantaneous blocking (IB) index:
#If LLB filter is true a third test criteria is applied that removes
#low latitude events that aren't true blocking because they don't block the flow.
#If MGI is True, return also the meridional gradient intensity, a metric of blocking strength
def IB_index(Z, lat_delta=15, lat0_min=30, lat0_max=75, thresh=-10, LLB_filter=False, LLB_thresh=-5, MGI=False):
    
    # set up memory storage
    IB = Z.sel(lat=slice(lat0_max, lat0_min)).copy()
    s_test = np.zeros([*IB.shape])
    n_test = np.zeros([*IB.shape])
    s2_test = np.zeros([*IB.shape])

    lons = IB.lon.values
    lats = IB.lat.values
    
    if MGI:
        mgi = np.zeros_like(s_test)
    
    # for each latitude and longitude check if conditions are true
    for j, lon in enumerate(lons):
        Zlon = Z.sel(lon=lon)
        for i, lat in enumerate(lats):
        
            # 1.) is the Southern gradient positive?
            ghgs = _GHGS(Zlon, lat+lat_delta, lat, lat-lat_delta, 0)
            s_test[:, i, j] = ghgs.values > 0
            
            # optionally return MGI
            if MGI:
                mgi[:, i, j] = ghgs.values
                
            # 2.) is the Northern gradient less than -10m/deg lat?
            ghgn = _GHGN(Zlon, lat+lat_delta, lat, lat-lat_delta, 0)
            n_test[:, i, j] = ghgn.values < thresh
            
            if LLB_filter:
                # 3.) is the Far southern gradient negative?
                ghgs2 = _GHGS2(Zlon, lat+lat_delta, lat, lat-lat_delta, 0)
                s2_test[:, i, j] = ghgs2.values < LLB_thresh
        
    # if both conditions are true than the test is passed:
    test = np.array(s_test) * np.array(n_test)
    
    if LLB_filter:
        test = test * np.array(s2_test)
        
    # pack the result in an xarray DataArray and return
    IB = xr.DataArray(test, coords=IB.coords, dims=IB.dims)
    IB.name = "IB index"
    IB.attrs['long_name'] = "IB index"
    
    if MGI:
        mgi[mgi < 0] = 0
        mgi_ix = IB.copy()
        mgi_ix.name = "meridional gradient index"
        mgi_ix.attrs['long_name'] = "meridional gradient index"
        mgi_ix.values = mgi
        return IB, mgi_ix
    
    return IB



#Compute Large Scale Blockings from instantaneous blockings
#based on whether they extend over at least +/-lon_thresh deg lat:
#(Assumes longitude is third coord)
def LSB_index(IB_ix, lon_thresh=7.5):
    
    lons = IB_ix.lon.values
    LSB_ix = IB_ix.copy()
    
    dL = lons[1] - lons[0]
    Nlats = np.ceil(lon_thresh / dL)
    
    # loop over longitudes
    for i, lon in enumerate(lons):
        LSBs = []
        # loop over offsets:
        for n in range(int(2*Nlats+1)):
            # extract the slice along longitude
            lon_slice = IB_ix.sel(lon=slice(lon-dL*n, lon+dL*((2*Nlats)-n)))
    
            # check if all points are blocked
            LSB = np.all(lon_slice.values, axis=2)
            LSBs.append(LSB)
        
        # if for any offset we find the whole slice contiguously blocked
        # then the point is part of a large scale block.
        LSB_ix.values[:, :, i] = np.any(np.array(LSBs), axis=0)
        
        
    LSB_ix.name = "LSB index"
    LSB_ix.attrs['long_name'] = "LSB index"
    return LSB_ix


# %%
