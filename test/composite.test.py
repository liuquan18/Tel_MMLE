#%%
import numpy as np
import pandas as pd
import xarray as xr
from src.composite.composite import reduce_var



#%%
def test_reduce_var():
    # create test data
    data = xr.DataArray(np.random.rand(10, 5, 5), dims=["com", "lat", "lon"])
    index = xr.DataArray(np.random.rand(10, 5, 5) > 0.5, dims=["com", "lat", "lon"])

    # test mean reduction
    mean_composite = reduce_var(index, data, dim="com", reduction="mean")
    assert mean_composite.shape == (5, 5)

    # test weighted mean reduction
    weights = xr.DataArray(np.random.rand(5), dims=["lat"])
    weighted_mean_composite = reduce_var(
        index, data, dim="com", reduction="mean_weighted", weights=weights
    )
    assert weighted_mean_composite.shape == (5, 5)

    # test mean reduction with bootstrap
    bootstrap_mean_composite = reduce_var(
        index, data, dim="com", reduction="mean", bootstrap=True
    )
    assert len(bootstrap_mean_composite) == 1000
    assert bootstrap_mean_composite[0].shape == (5, 5)

    # test mean reduction with same number of extremes
    mean_same_number_composite = reduce_var(
        index, data, dim="com", reduction="mean_same_number", count=3
    )
    assert mean_same_number_composite.shape == (5, 5)
# %%
