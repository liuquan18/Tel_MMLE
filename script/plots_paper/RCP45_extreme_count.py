# %%
import xarray as xr
import matplotlib.pyplot as plt
import src.plots.extreme_plot as extplt
import importlib

importlib.reload(extplt)
import proplot as pplt
import numpy as np
from scipy.stats import linregress

# %%


def read_extrc(model, fixed_pattern="decade_mpi"):
    """read extreme counts"""
    odir = "/work/mh0033/m300883/Tel_MMLE/data/" + model + "/extreme_count/"
    filename = f"plev_50000_{fixed_pattern}_first_JJA_extre_counts.nc"
    ds = xr.open_dataset(odir + filename).pc

    # divide the ensemble size of each model
    ens_sizes = {
        "MPI_GE": 100,
        "MPI_GE_onepct": 100,
        "MPI_GE_RCP45": 100,
        "CanESM2": 50,
        "CESM1_CAM5": 40,
        "MK36": 30,
        "GFDL_CM3": 20,
    }
    ds = ds / ens_sizes[model]
    return ds


def test_trend_significance(dataarray):
    time = np.arange(dataarray.time.size)
    values = dataarray.values
    slope, intercept, r_value, p_value, std_err = linregress(time, values)
    return p_value


def plot_extreme_counts(extrc, ax, extr_type, model_color, ci = True):
    
    line = extrc.sel(extr_type=extr_type,confidence = 'true').plot.line(
            ax=ax, 
            x = 'time',
            color = 'k',
            linewidth = 2,)

    if ci:
        # fill between the confidence interval ['low','high']
        ax.fill_between(
            extrc.time,
            extrc.sel(extr_type=extr_type,confidence = 'low').values,
            extrc.sel(extr_type=extr_type,confidence = 'high').values,
            color = 'gray',
            alpha = 0.3,
        )

    return line
# %%
extreme_counts = read_extrc("MPI_GE_RCP45")

# %%
count_pos = extreme_counts.sel(extr_type="pos", confidence="true", mode="NAO")
count_neg = extreme_counts.sel(extr_type="neg", confidence="true", mode="NAO")
# %%
count_pos_1950 = count_pos.sel(time=slice("1950", None))
count_neg_1950 = count_neg.sel(time=slice("1950", None))
# %%
# Test the significance of the trend in the two time series

p_value_pos_1950 = test_trend_significance(count_pos_1950)
p_value_neg_1950 = test_trend_significance(count_neg_1950)

# %%
p_value_pos = test_trend_significance(count_pos)
p_value_neg = test_trend_significance(count_neg)

# %%
extrcs = {"MPI_GE_RCP45": extreme_counts}
#%%
eof_rcp45 = xr.open_dataset("/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_RCP45/EOF_result/plev_50000_decade_mpi_first_JJA_eof_result.nc")
pc_first = eof_rcp45.pc.sel(time = slice('1850','1859'))
pc_last = eof_rcp45.pc.sel(time = slice('2090','2099'))

#%%
fig, axes = plt.subplots(1, 3, figsize=(12, 5))
pc_first.plot.hist(ax=axes[0], alpha=0.7, bins=np.arange(-4, 4.1, 0.5), color='k', label='first10')
pc_last.plot.hist(ax=axes[0], alpha=0.7, bins=np.arange(-4, 4.1, 0.5), color='r', label='last10')

# vline at mean
axes[0].axvline(pc_first.mean(), color='k', linestyle='dashed', linewidth=1.5)
axes[0].axvline(pc_last.mean(), color='r', linestyle='dashed', linewidth=1.5)


# Remove right and top spines
for ax in axes:
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

pos_line = plot_extreme_counts(extreme_counts.sel(mode = 'NAO'), axes[1], 'pos', 'k')
neg_line = plot_extreme_counts(extreme_counts.sel(mode = 'NAO'), axes[2], 'neg', 'k')

axes[1].set_ylabel('extreme counts')
axes[2].set_ylabel('extreme counts')

axes[0].set_title('')
axes[1].set_title('')
axes[2].set_title('')
axes[0].set_xlabel('NAO index')
# add a, b,c
axes[0].text(0.05, 0.95, 'a', transform=axes[0].transAxes, fontsize=16, fontweight='bold', va='top')
axes[1].text(0.05, 0.95, 'b', transform=axes[1].transAxes, fontsize=16, fontweight='bold', va='top')
axes[2].text(0.05, 0.95, 'c', transform=axes[2].transAxes, fontsize=16, fontweight='bold', va='top')


# add value of p-value as text
axes[1].text(
    0.95,
    0.95,
    f"1850-2090 p-value: {p_value_pos:.2f}",
    transform=axes[1].transAxes,
    fontsize=8,
    ha="right",
    va="top",
)

axes[1].text(
    0.95,
    0.9,
    f"1950-2090 p-value: {p_value_pos_1950:.2f}",
    transform=axes[1].transAxes,
    fontsize=8,
    ha="right",
    va="top",
)

axes[2].text(
    0.95,
    0.95,
    f"1850-2090 p-value: {p_value_neg:.2f}",
    transform=axes[2].transAxes,
    fontsize=8,
    ha="right",
    va="top",
)

axes[2].text(
    0.95,
    0.9,
    f"1950-2090 p-value: {p_value_neg_1950:.2f}",
    transform=axes[2].transAxes,
    fontsize=8,
    ha="right",
    va="top",
)
axes[1].set_ylabel('extreme counts')
axes[2].set_ylabel('extreme counts')

axes[0].set_title('')
axes[1].set_title('')
axes[2].set_title('')
axes[0].set_xlabel('NAO index')
axes[0].legend()
# all axes no grid
for ax in axes:
    ax.grid(False)
    ax.tick_params(axis='x', which='minor', bottom=False)
    ax.tick_params(axis='y', which='minor', left=False)
axes[0].set_xticks(np.arange(-3, 3.1, 1.5))



plt.savefig(
    "/work/mh0033/m300883/Tel_MMLE/docs/source/plots/paper_supplymentary/MPI_GE_RCP45_extreme_count.pdf",
    bbox_inches="tight",
    dpi=300,
)

# %%
