import numpy as np
import pandas as pd
import xarray as xr


def extreme(
    xarr: xr.DataArray,
    extreme_type: str,
    threshold: int = 2,
) -> xr.DataArray:
    """
    select only the positive or extreme cases from the index.
    detect the extreme cases identified from the threshold.
    **Arguments**
        *xarr* the index to be checked.
        *extrem_type*: 'pos' or 'neg
        *threshold* the threshold to identify extreme cases.
    **Return**
        *extreme* the extreme dataArray with neg and pos.
    """
    if extreme_type == "pos":
        extreme = xarr.where(xarr > threshold, drop=True)
    elif extreme_type == "neg":
        extreme = xarr.where(xarr < -1 * threshold, drop=True)

    return extreme


def composite(
    index: xr.DataArray, data: xr.DataArray, reduction: str = "mean", dim: str = "com"
):
    """
    do composite mean or count of data given the index
    **Argument**
        *index* the index which declears the coords to do the composite analysis
        *data* the data to be meaned or counted
        *reduction* "mean" or "count"
        *dim* along which dim to do the mean and count.
    """
    # get the data at the  coordinates
    sel_data = data.where(index)

    if reduction == "mean":
        composite = sel_data.mean(dim=dim)
    elif reduction == "count":
        composite = sel_data.count(dim=dim)
    return composite


def extreme_composite(index, data, reduction="mean", dim="com",threshold = 2):
    """
    the composite mean or count of data, in terms of different extreme type.

    **Arguments**
        *index* the from which the coordinates of
        extreme neg or pos cases are determined.
        *data* the field that are going to be selected and averaged.
    **Return**
        *extreme_composite* the mean field or counts of extreme cases.
    """
    Ext_composite = []
    extreme_type = xr.DataArray(["pos", "neg"], dims=["extr_type"])
    for extr_type in extreme_type.values:

        # get the coordinates of the extremes
        extr_index = extreme(index, extreme_type=extr_type,threshold=threshold)

        # do composite analysis based on the extreme index
        extr_composite = composite(extr_index, data, reduction=reduction, dim=dim)
        Ext_composite.append(extr_composite)

    # concate the positive and negative extremes together.
    Extreme_composite = xr.concat(Ext_composite, dim=extreme_type)

    return Extreme_composite
