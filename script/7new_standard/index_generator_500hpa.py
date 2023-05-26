#%%
#%%
import src.MMLE_TEL.index_generator as index_generate
# %%
MPIGE_decade = index_generate.decompose_plev("MPI_GE_onepct",plev = 50000,fixedPattern = 'decade',standard='temporal_ens')

# %%
MPIGE_decade.save_result()
# %%
MPIGE_all = index_generate.decompose_plev("MPI_GE_onepct",plev = 50000,fixedPattern = 'all',standard='temporal_ens')

# %%
MPIGE_all.save_result()
# %%
CESM_all = index_generate.decompose_plev("CESM1_CAM5",plev = 50000,fixedPattern = 'all',standard='temporal_ens')
# %%
CESM_all.save_result()

# %%
Can_all = index_generate.decompose_plev("CanESM2",plev = 50000,fixedPattern = 'all',standard='temporal_ens')
# %%
Can_all.save_result()
# %%

# %%
MPIGE_decade_trop = index_generate.decompose_troposphere("MPI_GE_onepct",vertical_eof='ind',fixedPattern = 'decade',standard='temporal_ens')
# %%
MPIGE_decade_trop.save_result()
# %%
# same for all
MPIGE_all_trop = index_generate.decompose_troposphere("MPI_GE_onepct",vertical_eof='ind',fixedPattern = 'all',standard='temporal_ens')
# %%
pc = MPIGE_all_trop.eof_result['pc']
std_pc = (pc - pc.mean(dim = ('time','ens')))/pc.std(dim = ('time','ens'))
MPIGE_all_trop.std_eof_result['pc'] = std_pc
MPIGE_all_trop.save_result()
# %%
