# %%
import xarray as xr
import numpy as np
from scipy.stats import skew
import matplotlib.pyplot as plt

# %%
# func to read eof data
def read_eof(model,fixedPattern):
    dir = f"/work/mh0033/m300883/Tel_MMLE/data/{model}/EOF_result/plev_50000_{fixedPattern}_temporal_ens_eof_result.nc"
    pc = xr.open_dataset(dir).pc
    return pc
# %%
def skewness(pc,mode):
    skews = []
    for i in range(0,pc.ens.size):
        skews.append(skew(pc.sel(ens = i,mode = mode)))
    return skews
# %%
# plot hist of the skewness
def plot_skewness(model,mode,fixedPattern = 'decade'):
    pc = read_eof(model,fixedPattern = 'decade')
    skews = skewness(pc,mode)
    plt.hist(skews,bins = np.arange(-0.5,0.51,0.1))
    plt.savefig(f"/work/mh0033/m300883/Tel_MMLE/docs/source/plots/skewness/{model}_skewness_{fixedPattern}.png")

# %%
# %%
plot_skewness('MPI_GE_onepct','NAO')
# %%
plot_skewness('MPI_GE','NAO')

# %%
plot_skewness('CanESM2','NAO')

# %%
plot_skewness('CESM1_CAM5','NAO')
# %%
# write the above plots (in'/work/mh0033/m300883/Tel_MMLE/docs/source/plots/skewness' ) into a md file
with open("/work/mh0033/m300883/Tel_MMLE/docs/source/skewness.md","w") as f:
    f.write("# Skewness of NAO\n")
    f.write("## MPI_GE_onepct\n")
    f.write("![MPI_GE_onepct](plots/skewness/MPI_GE_onepct_skewness_decade.png)\n")
    f.write("## MPI_GE\n")
    f.write("![MPI_GE](plots/skewness/MPI_GE_skewness_decade.png)\n")
    f.write("## CanESM2\n")
    f.write("![CanESM2](plots/skewness/CanESM2_skewness_decade.png)\n")
    f.write("## CESM1_CAM5\n")
    f.write("![CESM1_CAM5](plots/skewness/CESM1_CAM5_skewness_decade.png)\n")
# %%
CanESM2 = read_eof('CanESM2','decade')
# %%
CAN_skew = skewness(CanESM2,'NAO')
# %%
# get the index where the skewness is smaller than 0
Can_skew = np.array(CAN_skew)
indices = np.where(Can_skew < 0)
# %%
indices
# %%
