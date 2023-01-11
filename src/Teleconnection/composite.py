import xarray as xr
import pandas as pd
import numpy as np


def extreme(
    xarr: xr.DataArray,
    extreme_type: str,
    threshold: int=2,
)-> xr.DataArray:
    """
    select only the extreme cases from the index.
    detect the extreme cases identified from the threshold.
    **Arguments**
        *xarr* the index to be checked.
        *extrem_type*: 'pos' or 'neg
        *threshold* the threshold to identify extreme cases.
    **Return**
        *extreme* the extreme dataArray with neg and pos.
    """
    if extreme_type=='pos':
        extreme = xarr.where(xarr>threshold,drop=True)
    elif extreme_type == 'neg':
        extreme = xarr.where(xarr<-1*threshold,drop=True)
    
    return extreme

def extreme_index(
    index:xr.DataArray,
    threshod: int
)->xr.DataArray:
    """
    the pos and neg extreme index.
    given an index of NAO and EA, return the extreme index of pos and neg
    **Arguments**
        *index* the NAO and EA index
        *thredshold* the extreme threshold 
    **Return**
        *extr_index* the extreme pos and neg index
    """
    extr_index = []
    extreme_type = xr.DataArray(['pos','neg'],dims = ['extr_type'])
    for extr_type in extreme_type.values:
        exindex = extreme(index,extreme_type=extr_type,threshold=threshod)
        extr_index.append(exindex)
    extr_index = xr.concat(extr_index, dim=extreme_type)
    return extr_index


def _composite(index,data,reduction='mean'):
    """
    composite analysis for single height layer.
    the composite mean or count of data, determined by the extreme state
    of index.
    **Arguments**
        *index* the extreme index (single NAO or EA, pos or neg)
        *data* the field that are going to be selected and averaged.
    **Return**
        *extreme_composite* the mean field or counts of extreme cases.
    """
    if ('hlayers' in index.dims) | ('plev' in index.dims):
        try:
            data = data.sel(hlayers = index.hlayers)
        except KeyError:
            data = data

    data = data.stack(com = ('time','ens'))
    index = index.stack(com = ('time','ens'))
    extr_index = extreme_index(index)
    extr_data = data.where(extr_index)

    if reduction == 'mean':
        composite = extr_data.mean(dim = 'com')
    elif reduction == 'count':
        composite = extr_index.count(dim = 'com')

    return composite

def composite(
    index:xr.DataArray,
    data:xr.DataArray,
    dims: str=('mode','extr_type'),
    reduction: str = 'mean',
    threshod = 2):
    """
    composite analysis for mulity dim index
    **Arguments**
        *index* the NAO and EA index
        *data* the data to be compositely analized.
        *dims* which dims is reserved.
        *reduction* 'mean' or 'count'
    **return**
        *composite* the composite mean or count of data
    """
    extre_index = extreme_index(index,threshod)
    extre_index = extre_index.stack(com = dims)
    composite = extre_index.groupby('com').apply(_composite,
    data = data,reduction=reduction)
    return composite

def composite(
    index:xr.DataArray,
    data:xr.DataArray,
    reduction: str='mean',
    period: str= 'all'):
    """
    composite mean maps of extreme cases or counts of extreme cases.
    given the index and the data, determine the time and ens coordinates 
    where extreme cases happen from the index, and select the data with 
    those indexes, average the selected fields.
    **Arguments**
        *index* the index of NAO and EA
        *data* the original geopotential data.
        *reduction* mean or count
        *period* 'first10','last10','all'
    **Return**
        *compostie* the composite mean of the extreme cases.
    """
    if period == 'all':
        index = index
    elif period == 'first10':
        index = index.isel(time = slice(0,10)) # the index actually started from 1856,so wrong here.
    elif period == 'last10':
        index = index.isel(time = slice(-10,None))
    data  = data.sel(time = index.time)


    Composite = []
    for mode in index.mode:
        _index = index.sel(mode = mode)
        composite = _index.groupby('hlayers').apply(composite,data = data,reduction=reduction)
        Composite.append(composite)
    Composite = xr.concat(Composite,dim = index.mode)

    return Composite




    