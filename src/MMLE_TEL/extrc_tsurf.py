
# %%
def decadal_extrc_tsurf(index:xr.DataArray,temp:xr.DataArray,hlayers: int= 50000):
    """
    extract the extreme count and the mean surface temperature every ten years.
    **Arguments**
        *index* the DataArray of pc index
        *temp* the -fldmean, -yearmean -ensmean surface temperauter.
    **Return**
        *ext_counts* the DataArray of extreme count, *time, *extr_type, *mode
        *t_surf_mean* the mean t_surface,
        the time here use the first year of the decade.
    """

    index = index.sel(hlayers = hlayers)

    ext_counts = []
    t_surf_mean = []
    for i in range(0,index.time.size,10):

        period = slice(i,i+10)
    
        period_pc = index.isel(time = period)
        period_tm = temp.isel(time = period)

        # extreme count
        period_ext_count = extreme.period_extreme_count(period_pc)
        period_mean_t = period_tm.mean(dim = 'time')

        ext_counts.append(period_ext_count)
        t_surf_mean.append(period_mean_t)

    periods = index.isel(time = np.arange(0,index.time.size,10)).time
    ext_counts = xr.concat(ext_counts, dim=periods)
    t_surf_mean = xr.concat(t_surf_mean, dim = periods)
    return ext_counts,t_surf_mean