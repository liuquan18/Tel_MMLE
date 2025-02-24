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
# %%
fig, axes = pplt.subplots(
    nrows=1,
    ncols=2,
    width=150 / 25.4,
    height=100 / 25.4,
)
fig.format(
    abc=True,
)

ax1 = axes[0]
ax2 = axes[1]

ax1, line_pos = extplt.extrc_time_line_single(
    extrcs,
    extr_type="pos",
    ax=ax1,
    ylim=(1.5, 4.5),
    ci=True,
    models=["MPI_GE_RCP45"],
)

ax2, line_neg = extplt.extrc_time_line_single(
    extrcs,
    extr_type="neg",
    ax=ax2,
    ylim=(1.5, 4.5),
    ci=True,
    models=["MPI_GE_RCP45"],
)
# add value of p-value as text
ax1.text(
    0.95,
    0.95,
    f"1850-2090 p-value: {p_value_pos:.2f}",
    transform=ax1.transAxes,
    fontsize=8,
    ha="right",
    va="top",
)

ax1.text(
    0.95,
    0.9,
    f"1950-2090 p-value: {p_value_pos_1950:.2f}",
    transform=ax1.transAxes,
    fontsize=8,
    ha="right",
    va="top",
)

ax2.text(
    0.95,
    0.95,
    f"1850-2090 p-value: {p_value_neg:.2f}",
    transform=ax2.transAxes,
    fontsize=8,
    ha="right",
    va="top",
)

ax2.text(
    0.95,
    0.9,
    f"1950-2090 p-value: {p_value_neg_1950:.2f}",
    transform=ax2.transAxes,
    fontsize=8,
    ha="right",
    va="top",
)


ax1.format(
    xlabel="",
    xlabelpad=0.8,
    xtickminor=False,
    xrotation=45,
    ylabel="Extreme occurence decade$^{-1}$ realization$^{-1}$",
    ylim=(1.5, 3.6),
    facecolor="none",
)
ax1.set_yticks(np.arange(1.5, 3.5, 0.5))

ax1.spines["right"].set_visible(False)
ax1.spines["top"].set_visible(False)

#### ax4 ####
# set the axis
ax2.spines["right"].set_visible(False)
ax2.spines["top"].set_visible(False)
ax2.format(
    xlabel="",
    xlabelpad=0.8,
    xtickminor=False,
    ylabel="",
    ylim=(1.5, 3.6),
    xrotation=45,
    facecolor="none",
)

plt.savefig(
    "/work/mh0033/m300883/Tel_MMLE/docs/source/plots/paper_supplymentary/MPI_GE_RCP45_extreme_count.pdf",
    bbox_inches="tight",
    dpi=300,
)

# %%
