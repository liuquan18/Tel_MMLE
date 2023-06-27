#%%
import xarray as xr
import src.MMLE_TEL.story_line as story_line
import src.plots.extrc_tsurf_scatter as extrc_tsurf
import src.MMLE_TEL.index_stats as index_stats
import src.plots.extrc_tsurf_scatter as extrc_tsurf

# %%
MJJA_first = xr.open_dataset("/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_onepct/EOF_result/troposphere_ind_decade_first_MJJA_eof_result.nc")

#%%
# read tsurf
def read_tsurf(tsurf_dir):
    tsurf = xr.open_dataset(tsurf_dir)

    # read the array either tusrf, ts, or tas
    try:
        tsurf_arr = tsurf.tsurf.squeeze()
    except AttributeError:
        try:
            tsurf_arr = tsurf.ts.squeeze()
        except AttributeError:
            tsurf_arr = tsurf.tas.squeeze()
    # change the temp time into datetime64
    try:
        tsurf_arr["time"] = tsurf_arr.indexes["time"].to_datetimeindex()
    except AttributeError:
        pass
    tsurf_arr = tsurf_arr - tsurf_arr[0]
    return tsurf_arr
# %%
temperature = read_tsurf("/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_onepct/ts_processed/ens_fld_year_mean.nc")
# %%
pcs = MJJA_first.pc.sel(plev = 85000)
# %%
