# %%
import src.Teleconnection.spatial_pattern as ssp
# %%

# %%
def decompose_period(xarr, nmode = 2,period = True):
    xarr = xarr.fillna(0) # fill nan with 0
    field = xarr.sortby("time")
    if 'ens' in xarr.dims:
        field = field.stack(com=("ens", "time"))
        dim = 'com'
    else:
        dim = 'time'
    standard = 'eof_spatial_std' if period else 'pc_temporal_std'
    eof_result = ssp.doeof(field, standard=standard,nmode=nmode,dim = dim)
    return eof_result


# %%
def standard_period(first_eof, last_eof):
    if 'ens' in first_eof.dims:
        mean = first_eof["pc"].mean(dim=("time", "ens"))
        std = first_eof["pc"].std(dim=("time", "ens"))
    else:
        mean = first_eof["pc"].mean(dim="time")
        std = first_eof["pc"].std(dim="time")

    first_eof["pc"] = (first_eof["pc"] - mean) / std
    last_eof["pc"] = (last_eof["pc"] - mean) / std
    return first_eof, last_eof