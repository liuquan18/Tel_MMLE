# %%
import xarray as xr
import numpy as np
from scipy import stats
# %%
############### Read EOFs ###############
def read_eof_decade(model, fixed_pattern="decade_mpi"):
    """read eofs that is decomposed by decade"""
    odir = f"/work/mh0033/m300883/Tel_MMLE/data/{model}/EOF_result/"
    filename = f"plev_50000_{fixed_pattern}_first_JJA_eof_result.nc"
    ds = xr.open_dataset(odir + filename)
    ds = ds.sel(mode="NAO")
    return ds

def split_first_last(eof_result):
    times = eof_result.time
    years = np.unique(times.dt.year)
    first_years = years[:10]
    last_years = years[-10:]

    eof_first = eof_result.isel(decade=0).sel(
        time=eof_result["time.year"].isin(first_years)
    )
    eof_last = eof_result.isel(decade=-1).sel(
        time=eof_result["time.year"].isin(last_years)
    )
    return eof_first.pc.values.flatten(), eof_last.pc.values.flatten()

def read_eof_rean(model, group_size=40):
    odir = "/work/mh0033/m300883/Tel_MMLE/data/" + model + "/"
    first_eof_path = odir + "EOF_result/first_plev50000_eof.nc"
    last_eof_path = odir + "EOF_result/last_plev50000_eof.nc"

    first_pc = xr.open_dataset(first_eof_path).pc.values.flatten()
    last_pc = xr.open_dataset(last_eof_path).pc.values.flatten()
    return first_pc, last_pc

#%%
MPI_GE_first_pc, MPI_GE_last_pc = split_first_last(read_eof_decade("MPI_GE"))
CanESM2_first_pc, CanESM2_last_pc = split_first_last(read_eof_decade("CanESM2"))
CESM1_CAM5_first_pc, CESM1_CAM5_last_pc = split_first_last(read_eof_decade("CESM1_CAM5"))
MK36_first_pc, MK36_last_pc = split_first_last(read_eof_decade("MK36"))
GFDL_CM3_first_pc, GFDL_CM3_last_pc = split_first_last(read_eof_decade("GFDL_CM3"))

CR20_first_pc, CR20_last_pc = read_eof_rean("CR20_allens")

############ KS test ############
# %%
def KS_test(first_pc, last_pc, whole = True):
    ks_res = stats.kstest(first_pc, last_pc, alternative='two-sided')
    return ks_res.pvalue


# %%
MPI_p = KS_test(MPI_GE_first_pc, MPI_GE_last_pc)
CanESM2_p = KS_test(CanESM2_first_pc, CanESM2_last_pc)
CESM1_CAM5_p = KS_test(CESM1_CAM5_first_pc, CESM1_CAM5_last_pc)
MK36_p = KS_test(MK36_first_pc, MK36_last_pc)
GFDL_CM3_p = KS_test(GFDL_CM3_first_pc, GFDL_CM3_last_pc)
CR20_p = KS_test(CR20_first_pc, CR20_last_pc)


# %%
MPI_p, CanESM2_p, CESM1_CAM5_p, MK36_p, GFDL_CM3_p, CR20_p



# %%
################# bootstrap test std #################
def bootstrap(first_pc, last_pc, n_resamples=1000, alpha=0.05):
    # calculate the difference in the std between first_pc and last_pc, and repeat for n_resamples
    # calculate the quantile between alpha/2 and 1-alpha/2
    # if the difference in the std is outside the quantile, then the difference is significant
    std_diffs = []
    for _ in range(n_resamples):
        resample_first = np.random.choice(first_pc, size=len(first_pc), replace=True)
        resample_last = np.random.choice(last_pc, size=len(last_pc), replace=True)
        std_diff = np.std(resample_first) - np.std(resample_last)
        std_diffs.append(std_diff)

    low_bnd = alpha / 2.0
    high_bnd = 1 - alpha / 2.0

    ci = np.quantile(std_diffs, [low_bnd, high_bnd])

    # check if 0 is in the interval, return boolean
    diff_significant = (ci[0] > 0) or (ci[1] < 0)

    return diff_significant

# %%
MPI_bootstrap = bootstrap(MPI_GE_first_pc, MPI_GE_last_pc)
CanESM2_bootstrap = bootstrap(CanESM2_first_pc, CanESM2_last_pc)
CESM1_CAM5_bootstrap = bootstrap(CESM1_CAM5_first_pc, CESM1_CAM5_last_pc)
MK36_bootstrap = bootstrap(MK36_first_pc, MK36_last_pc)
GFDL_CM3_bootstrap = bootstrap(GFDL_CM3_first_pc, GFDL_CM3_last_pc)
CR20_bootstrap = bootstrap(CR20_first_pc, CR20_last_pc)
# %%
MPI_bootstrap, CanESM2_bootstrap, CESM1_CAM5_bootstrap, MK36_bootstrap, GFDL_CM3_bootstrap, CR20_bootstrap
# %%
# %%
