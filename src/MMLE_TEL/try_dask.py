import dask
import dask.distributed
import xarray as xr
import src.extreme.extreme_ci as extreme

# Open the input data files as xarray.Dataset objects
eof_result = xr.open_dataset("/work/mh0033/m300883/Tel_MMLE/data/CanESM2/EOF_result/plev_50000_decade_temporal_ens_eof_result.nc").pc
tsurf_mean = xr.open_dataset("/work/mh0033/m300883/Tel_MMLE/data/CanESM2/ts_processed/ens_fld_year_mean.nc").ts.squeeze()
tsurf_increase = tsurf_mean - tsurf_mean[0]

# Define the index and tsurf variables
index = eof_result
tsurf = tsurf_increase

# Create a Dask client to distribute the computation
client = dask.distributed.Client()

ext_counts = []
t_surf_mean = []

# Define a function to compute the extreme count and mean temperature for a single period
@dask.delayed
def compute_period(period):
    period_pc = index.isel(time=period)
    period_tm = tsurf.sel(time=period_pc.time, method='nearest')
    time_tag = period_pc.time[0] # for reference 
    print(time_tag.dt.year.values)

    # extreme count
    period_ext_count = extreme.extreme_count_xr(period_pc, ci=ci)
    period_mean_t = period_tm.mean(dim="time")

    period_ext_count['time'] = time_tag
    period_mean_t['time'] = time_tag

    # set time as the new dimension
    period_ext_count = period_ext_count.expand_dims('time')
    period_mean_t = period_mean_t.expand_dims('time')

    return period_ext_count, period_mean_t

# Compute the extreme count and mean temperature for each period in parallel
results = []
for i in range(0, index.time.size, 10): # avoid the last period which is not 10 years
        
    period = slice(i, i + 10)
    result = compute_period(period)
    results.append(result)

ext_counts, t_surf_mean = dask.compute(*results)

# Close the Dask client
client.close()