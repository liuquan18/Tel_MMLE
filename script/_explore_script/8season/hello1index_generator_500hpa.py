# %%
# %%
import src.MMLE_TEL.index_generator as index_generate
import xarray as xr
import numpy as np
import cartopy.crs as ccrs
import src.compute.slurm_cluster as scluster
import concurrent.futures
import sys
from mpi4py import MPI

# %%
import importlib

importlib.reload(index_generate)
importlib.reload(scluster)


# %%
# function for generate the index
def index_gen(model, fixedPattern = 'decade', plev=50000, season="MJJA", standard="first"):
    generator = index_generate.decompose_plev(
        model, plev=plev, fixedPattern=fixedPattern, standard=standard, season=season
    )
    generator.save_result()


def main():
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    models = {
        0: "MPI_GE",
        1: "CanESM2",
        2: "CESM1_CAM5",
        3: "MK36",
        4: "GFDL_CM3",
        # 5: "MPI_GE_onepct"
    }

    seasons = {
        0: "DJFM",
        1: "JJAS",
        2: "DJF",
        3: "MAM",
        4: "JJA",
        5: "SON",
    }

    print(f"rank {rank} of {size}")
    modelidx = rank // 6
    seasonidx = rank % 6

    model = models[modelidx]
    season = seasons[seasonidx]
    index_gen(model, season=season)
    print(f"rank {rank} finished")
    comm.Barrier()
#%%
    
if __name__ == "__main__":
    main()



"""

# %%


# for different models
# %%
def index_gen_models(season="MAM"):
    with concurrent.futures.ProcessPoolExecutor() as executor:
        futures = [
            # executor.submit(index_gen, 'MPI_GE_onepct', 'decade', plev=50000),
            executor.submit(
                index_gen,
                "MPI_GE",
                "decade",
                plev=50000,
                season=season,
                standard="first",
            ),
            executor.submit(
                index_gen,
                "CanESM2",
                "decade",
                plev=50000,
                season=season,
                standard="first",
            ),
            executor.submit(
                index_gen,
                "CESM1_CAM5",
                "decade",
                plev=50000,
                season=season,
                standard="first",
            ),
            executor.submit(
                index_gen, "MK36", "decade", plev=50000, season=season, standard="first"
            ),
            # executor.submit(index_gen, 'GFDL_CM3', 'decade', plev=50000,season = season,standard = 'first'),
        ]
        for future in concurrent.futures.as_completed(futures):
            try:
                result = future.result()
            except Exception as e:
                print(f"Exception: {e}")


# for different seasons and standards, only for MPI_GE_onepct
# %%
def index_gen_seasons(season="MAM"):
    with concurrent.futures.ProcessPoolExecutor() as executor:
        futures = [
            executor.submit(
                index_gen,
                "MPI_GE_onepct",
                "decade",
                plev=50000,
                season=season,
                standard="first",
            ),
            executor.submit(
                index_gen,
                "MPI_GE_onepct",
                "decade",
                plev=50000,
                season=season,
                standard="first",
            ),
            # executor.submit(index_gen, 'MPI_GE_onepct', 'decade', plev=50000, season='DJFM', standard='first'),
            # executor.submit(index_gen, 'MPI_GE_onepct', 'decade', plev=50000, season='JJAS', standard='temporal_ens'),
            # # executor.submit(index_gen, 'MPI_GE_onepct', 'decade', plev=50000, season='DJFM', standard='temporal_ens')
            # executor.submit(index_gen, 'MPI_GE_onepct', 'decade', plev=50000, season='MAM', standard='temporal_ens')
        ]
        for future in concurrent.futures.as_completed(futures):
            try:
                result = future.result()
            except Exception as e:
                print(f"Exception: {e}")


# define main funtion, run index_gen_seasons() or index_gen_models() here
# %%
def main():
    index_gen_models(season="MAM")
    # index_gen_models()


# run mean function here
if __name__ == "__main__":
    main()

"""