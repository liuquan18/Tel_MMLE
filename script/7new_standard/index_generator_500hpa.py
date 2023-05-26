#%%
#%%
import src.MMLE_TEL.index_generator as index_generate
import xarray as xr
# %%
def standard(generator):
    pc = generator.eof_result['pc']
    std_pc = (pc - pc.mean(dim = ('time','ens')))/pc.std(dim = ('time','ens'))
    generator.std_eof_result['pc'] = std_pc
    return generator

# %%
# function for generate the index
def index_gen(model,fixedPattern):
    generator = index_generate.decompose_plev(model,plev = 50000,fixedPattern = fixedPattern,standard='temporal_ens')
    decade = standard(generator)
    decade.save_result()
# %%
# CanESM2
index_gen('CanESM2','decade')
index_gen('CanESM2','all')
# %%
# CESM1_CAM5
index_gen('CESM1_CAM5','decade')
index_gen('CESM1_CAM5','all')
# %%
# MK3.6
index_gen('MK36','decade')
index_gen('MK36','all')
# %%
def read_eof(model,fixedPattern):
    eof_dir = f"/work/mh0033/m300883/Tel_MMLE/data/{model}/EOF_result/plev_50000_{fixedPattern}_temporal_ens_eof_result.nc"
    eof_result = xr.open_dataset(eof_dir)
    return eof_result
# %%
Can_decade_eof = read_eof('CanESM2','decade')
Can_all_eof = read_eof('CanESM2','all')
# %%
# CESM1_CAM5
CESM_decade_eof = read_eof('CESM1_CAM5','decade')
CESM_all_eof = read_eof('CESM1_CAM5','all')
# %%
# MK3.6
MK36_decade_eof = read_eof('MK36','decade')
MK36_all_eof = read_eof('MK36','all')
# %%
