#%%
import xarray as xr
import numpy as np
from src.blocking.utils import reg_lens
# %%

#An auxilliary function used by blocking_event_index
#Takes an index, and checks if any points in a lat lon box of +/-lat_thresh
# +/- lon thresh around each point is true. If so make that point true.
def _box_ix(ix, lat_inter = 5, lon_thresh=10):

    # rolling the index over a lat lon box of +/-lat_thresh, and construct
    # data array with the window dimension as lat_zone and lon_zone
    ix_r = ix.rolling(
        lat = lat_inter,
        lon = lon_thresh,
        center = True
    ).construct(window_dim = {"lat":"lat_zone", "lon":"lon_zone"})

    # fill the nan values at the edges with 0
    ix_r = ix_r.fillna(0)

    # assign the lat_zone and lon_zone coordinates
    ix_r = ix_r.assign_coords(
        lat_zone = np.arange(lat_inter) + 1,
        lon_zone = np.arange(lon_thresh) + 1)
    
    # check if any value in the box is true
    ix_box = ix_r.any(dim = ["lat_zone", "lon_zone"]).astype(np.int32)

    return ix_box


#An auxilliary function used by blocking_event_index.
# Takes in an index and only returns true when the index
#is true for pers_thresh consecutive timesteps or more.
#Assumes TxLatxLon cube


#%%
def _keep_point(ix, pers_thresh=5):
        L, S = reg_lens (ix)    
        keepS = S == 1
        keepL = L > pers_thresh
        keep_points = (np.repeat(keepL, L) * np.repeat(keepS, L))
        return keep_points 

# apply the _keep_point function to the whole index along time dimension
def _check_persistence(ix, pers_thresh=5):
    pers_ix = xr.apply_ufunc(
        _keep_point, 
        ix, 
        input_core_dims=[["time"]], 
        output_core_dims=[["time"]], 
        vectorize=True, 
        dask="parallelized", 
        output_dtypes=[np.int32],
        kwargs={"pers_thresh": pers_thresh}
    )
    return pers_ix

def _yearly_persistence(ix, pers_thresh=5):
    year_per = ix.groupby("time.year").apply(_check_persistence, pers_thresh=pers_thresh)
    return year_per


# %%
#Defines a blocking event by requiring spatial and temporal persistence:
def blocking_event_index(LSB_ix,lat_thresh=2.5,lon_thresh=5,pers_thresh=5):
    
    #Find all points where there is LSB within
    #a certain lat lon region
    box_index=_box_ix(LSB_ix,lat_thresh,lon_thresh)
    
    blocking_event=_yearly_persistence(box_index,pers_thresh)
    
    return blocking_event

#%%
