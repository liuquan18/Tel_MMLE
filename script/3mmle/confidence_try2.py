# do eof anaylsis on the xarray dataarray with dimension (time,lat,lon)
# using eof package
def doeof(xarr):
    """
    xarr: xarray dataarray with dimension (time,lat,lon)
    return: xarray dataarray with dimension (mode,lat,lon)
    """
    # convert the xarray dataarray to numpy array
    # with dimension (time,lat,lon)
    data = xarr.values
    # create a mask to ignore the nan values
    mask = np.isnan(data)
    # create a solver object
    solver = Eof(data,mask=mask)
    # calculate the principal components
    pc = solver.pcs(npcs=3,pcscaling=1)
    # calculate the eigenvalues
    eig = solver.eigenvalues()
    # calculate the variance explained by each mode
    var = solver.varianceFraction()
    # calculate the spatial pattern of each mode
    # with dimension (mode,lat,lon)
    pattern = solver.eofs(neofs=3)
    # create a xarray dataarray with dimension (mode,lat,lon)
    pc_xarr = xr.DataArray(pc,dims=['mode','lat','lon'],
                           coords={'mode':range(1,4),
                                   'lat':xarr.lat,
                                   'lon':xarr.lon})
    eig_xarr = xr.DataArray(eig,dims=['mode'],
                            coords={'mode':range(1,4)})
    var_xarr = xr.DataArray(var,dims=['mode'],
                            coords={'mode':range(1,4)})
    pattern_xarr = xr.DataArray(pattern,dims=['mode','lat','lon'],
                                coords={'mode':range(1,4),
                                        'lat':xarr.lat,
                                        'lon':xarr.lon})
    return pc_xarr,eig_xarr,var_xarr,pattern_xarr