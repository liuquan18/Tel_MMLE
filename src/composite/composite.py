import numpy as np
import pandas as pd
import xarray as xr


def extreme(
    xarr: xr.DataArray,
    extreme_type: str,
    threshold: float = 1.5,
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


def reduce_var(
    index: xr.DataArray,
    data: xr.DataArray,
    dim: str = "com",
    reduction="mean",
    bootstrap=False,
    **kwargs
):
    """select the var data based on the coordinates of extreme cases in index, and do the average along dim.
    if bootstrap is True, then do the bootstrap resampling.
    **Arguments**
        *index* the index of NAO and EA extremes
        *data* the ts or pr data to do average.
        *dim* the dimension to do average.
        *reduction* the method to do average, 'mean' or 'mean_same_number', 'mean_weighted'
        *bootstrap* if True, do the bootstrap resampling.
        *kwargs* the arguments for different reduction methods.for example, 'count' for 'mean_same_number'
    **Return**
        *composite_mean* the composite mean of the extreme cases.
    """
    if reduction == "mean":
        # get the data at the  coordinates
        sel_data = data.where(index)
        composite_mean = sel_data.mean(dim=dim)

    elif reduction == "mean_same_number":
        # get the number of extremes to reduce from kwargs
        num = kwargs.get("count", 'all')
        # sel the largest 40 values from index
        index = index.copy()
        index = index.squeeze()
        if index.attrs["extreme_type"] == "pos":
            index = index.sortby(index, ascending=False)
        elif index.attrs["extreme_type"] == "neg":
            index = index.sortby(index, ascending=True)
        index = index.isel(com=slice(0, num))  # select the first {num} indexes
        sel_data = data.where(index)
        composite_mean = sel_data.mean(dim=dim)


    elif reduction == "mean_weighted":
        sel_data = data.where(index)
        weights = index
        composite_mean = sel_data.weighted(weights).mean(dim=dim)
    
    # set count info as attribute
    composite_mean.attrs["reduction"] = reduction
    composite_mean.attrs["count"] = num

    if bootstrap:
        n_resamples = kwargs.get("n_resamples", 1000)
        # get the 'com' random index with the shape sel_data.size['com'] and 1000 times
        n_samples = sel_data.sizes[
            dim
        ]  # the resampled data is the same length as the original data
        rng = np.random.default_rng(seed=12345)
        sampled_index = rng.choice(n_samples, size=(n_samples, n_resamples), replace=True)
        composite_mean = []
        for i in range(1000):
            sample = sel_data.isel(com=sampled_index[:, i])
            composite_mean.append(sample.mean(dim=dim))

        composite_mean = xr.concat(composite_mean, dim="bootstrap")

    return composite_mean


def extreme_composite(
    index, data, reduction="mean", dim="com", threshold=1.5, bootstrap=False, **kwargs
):
    """
    get the indexes for the extreme cases, and do the composite analysis.
    **Arguments**
        *index* the index of NAO and EA 
        *data* the data to do average
        *reduction* the method to do average, 'mean' or 'mean_same_number', 'mean_weighted'
        *dim* the dimension to do average.
        *threshold* the threshold to identify extreme cases.
        *bootstrap* if True, do the bootstrap resampling.
        *kwargs* the arguments for different reduction methods.for example, 'count' for 'mean_same_number'
    **Return**
        *composite_mean* the composite mean of the extreme positive and negative cases.
    """
    Ext_composite = []
    extreme_type = xr.DataArray(["pos", "neg"], dims=["extr_type"])
    for extr_type in extreme_type.values:
        # get the coordinates of the extremes
        extr_index = extreme(index, extreme_type=extr_type, threshold=threshold)
        extr_index.attrs["extreme_type"] = extr_type

        # do composite analysis based on the extreme index
        extr_composite = reduce_var(
            extr_index,
            data,
            reduction=reduction,
            dim=dim,
            bootstrap=bootstrap,
            **kwargs
        )
        Ext_composite.append(extr_composite)

    # concate the positive and negative extremes together.
    Extreme_composite = xr.concat(Ext_composite, dim=extreme_type)

    return Extreme_composite


def Tel_field_composite(
    index: xr.DataArray,
    data: xr.DataArray,
    # threshold: float = 1.5,
    # reduction="mean",
    threshold,
    reduction,
    bootstrap=False,
    count = None,
):
    """
    composite mean maps of NAO and EA extremes.
    **Arguments**
        *index* the index of NAO and EA
        *data* the original geopotential data.
        *reduction* mean or count
        *period* 'first10','last10','all'
    **Return**
        *compostie* the composite mean of the extreme cases of NAO and EA.
    """

    # Select the same time period
    index_c = index.copy()  # make a copy of the original data
    data_c = data.copy()

    # select the same time period
    data_c = data_c.sortby("time")
    data_c = data_c.sel(time=index_c.time, method="nearest")
    data_c["time"] = index_c.time  # make sure the time is the same to use where function

    # combine time and ens into one dim
    index_c = index_c.stack(com=("time", "ens"))
    data_c = data_c.stack(com=("time", "ens"))

    # since there is 'mode' dim in index, here groupby.
    tel_composite = index_c.groupby("mode").apply(
        extreme_composite,
        data=data_c,
        reduction=reduction,
        dim="com",
        threshold=threshold,
        bootstrap=bootstrap,
        count = count,
    )

    return tel_composite


def first_last_extreme_composite(
    first_index: xr.DataArray,
    last_index: xr.DataArray,
    var_data: xr.DataArray,
    threshold: float = 1.5,
    reduction: str = "mean",
    return_diff: bool = True,
    alpha = 0.05,
    **kwargs
):
    """
    composite mean maps of first and last 10 years of field in terms of teleconnection mode extremes.
    **Arguments**
        *first_index* the index of first 10 years of NAO and EA
        *last_index* the index of last 10 years of NAO and EA
        *var_data* the ts or pr data.
        *reduction* the method to do the reduction, 'mean' or 'mean_same_number', 'mean_weighted'
        *return_diff* if True, return the difference between the first and last composite, and the significance.
    **Return**
        *composite* the composite mean of the extreme cases of NAO and EA for the first and last 10 years, and the difference between them.
    """
    count = kwargs.get('count', None)    # Select the same time period
    try:
        first_index = first_index.drop_vars(('decade','plev'))
        last_index = last_index.drop_vars(('decade','plev'))
    except ValueError:
        pass

    # first 10 years
    first_composite = Tel_field_composite(
        first_index,
        var_data,
        threshold=1.5,
        reduction=reduction,
        bootstrap=False,
        count = count,
    )

    # last 10 years
    last_composite = Tel_field_composite(
        last_index,
        var_data,
        threshold=threshold,
        reduction=reduction,
        bootstrap=False,
        count = count,
    )

    # combine the first and last composite
    period = xr.IndexVariable('period', ['first', 'last'])
    composite = xr.concat([first_composite, last_composite], dim=period)
    diff = last_composite - first_composite

    if return_diff:
        print(" doing bootstrap resampling...")
        # first 10 years with bootstrap
        first_composite_boot = Tel_field_composite(
            first_index,
            var_data,
            threshold=1.5,
            reduction=reduction,
            bootstrap=True,
        )

        # last 10 years
        last_composite_boot = Tel_field_composite(
            last_index,
            var_data,
            threshold=threshold,
            reduction=reduction,
            bootstrap=True,
        )

        # difference between first and last 10 years
        diff_boot = last_composite_boot - first_composite_boot

        # check if the difference is significant
        # get the 95% confidence interval
        low_bnd = alpha/2.0
        high_bnd = 1-alpha/2.0
        ci = diff_boot.quantile([low_bnd, high_bnd], dim="bootstrap")
        # check if 0 is in the interval, return boolean
        diff_sig = xr.where((ci.sel(quantile=low_bnd) > 0) | (ci.sel(quantile=high_bnd) < 0), 1, 0)
        # combine the first, last, diff and diff_sig together
        period = xr.IndexVariable('period', ['first', 'last', 'diff', 'diff_sig'])
        composite = xr.concat([first_composite, last_composite, diff, diff_sig], dim=period)

    return composite
