# %%
import src.MMLE_TEL.story_line as story_line
import numpy as np
import xarray as xr
import importlib
import matplotlib.pyplot as plt
import src.plots.extreme_plot as extrc_tsurf

# %%
# %%
def read_data( ens_size,model = 'MPI_GE'):
    eof_dir = f"/work/mh0033/m300883/Tel_MMLE/data/{model}_random/EOF_result/plev_50000_decade_temporal_ens{str(ens_size)}_eof_result.nc"
    eof_result = xr.open_dataset(eof_dir)

    tsurf_dir = f"/work/mh0033/m300883/Tel_MMLE/data/{model}_random/ts_processed/ens_fld_year_mean.nc"
    tsurf = xr.open_dataset(tsurf_dir).tsurf

    return eof_result, tsurf
# %%
# extreme event count vs. tsurf
def random_scatter(ci = 'bootstrap',ens_size = 20,model = 'MPI_GE'):
    eof_result,tsurf = read_data(ens_size,model = model)
    print("ploting the extreme event count vs. tsurf")
    try:
        tsurf_mean = tsurf.mean(dim="ens").squeeze()
    except ValueError:
        tsurf_mean = tsurf
    tsurf_increase = tsurf_mean - tsurf_mean[0]

    ext_counts, t_surf_mean = extrc_tsurf.decadal_extrc_tsurf(
        eof_result.pc, tsurf_increase,ci = ci
    )
    extrc_tsurf_scatter = extrc_tsurf.extCount_tsurf_scatter(
        ext_counts, t_surf_mean, 
    )
    plt.savefig(f"/work/mh0033/m300883/Tel_MMLE/docs/source/plots/random_ens/{model}_ens_{str(ens_size)}" + f"_{ci}_extreme_count_tsurf.png", dpi=300)

# %%
import concurrent.futures

with concurrent.futures.ProcessPoolExecutor() as executor:
    futures = []
    futures.append(executor.submit(random_scatter, ens_size=20,model = 'MPI_GE_onepct'))
    futures.append(executor.submit(random_scatter, ens_size=30,model = 'MPI_GE_onepct'))
    futures.append(executor.submit(random_scatter, ens_size=40,model = 'MPI_GE_onepct'))
    futures.append(executor.submit(random_scatter, ens_size=50,model = 'MPI_GE_onepct'))

    for future in concurrent.futures.as_completed(futures):
        try:
            result = future.result()
        except Exception as exc:
            print(f'generated an exception: {exc}')

# %%
