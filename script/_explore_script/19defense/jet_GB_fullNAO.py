# %%
import xarray as xr
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import os
from src.mechnisms.mechisms import read_greenland_blocking, read_jetStream, Jet_location


# %%
model = "MPI_GE"
ppc_first_ano = xr.open_dataset(
    f"/work/mh0033/m300883/Tel_MMLE/data/{model}/full_signal/first_pattern_projected_ano.nc"
).pc
ppc_last_ano = xr.open_dataset(
    f"/work/mh0033/m300883/Tel_MMLE/data/{model}/full_signal/last_pattern_projected_ano.nc"
).pc
# %%

ppc_first_ano = ppc_first_ano.isel(mode=0)
ppc_last_ano = ppc_last_ano.isel(mode=0)


# %%
def anomaly(x):
    return x - x.sel(time=slice("1850", "1859")).mean()


# %%
def corresponding_NAO(NAO, var, threshold=0):
    """
    var: jet_stream_loc or GB
    """
    std = var.sel(time=slice("1850", "1859")).std()
    var_pos = var.where(var > threshold * std)
    NAO_var_pos = NAO.where(var_pos.notnull())

    return NAO_var_pos


# %%
# length in cm
cm = 1 / 2.54
fig, axes = plt.subplots(1, 1, figsize=(34 * cm, 19 * cm))

ppc_first_ano.plot.hist(
    ax=axes,
    alpha=0.7,
    bins=np.arange(-4, 4.1, 0.5),
    color="k",
    label="1850-1859",
    edgecolor="white",
    linewidth=0.5,
)
ppc_last_ano.plot.hist(
    ax=axes,
    alpha=0.7,
    bins=np.arange(-4, 4.1, 0.5),
    color="r",
    label="2090-2099",
    edgecolor="white",
    linewidth=0.5,
)

# vline at mean
axes.axvline(ppc_first_ano.mean(), color="k", linestyle="dashed", linewidth=3)
axes.axvline(ppc_last_ano.mean(), color="r", linestyle="dashed", linewidth=3)

# fontsize for ticks and labels
axes.tick_params(axis="both", which="major", labelsize=14)
axes.set_xlabel("standarized NAO index", fontsize=20)
axes.set_ylabel("Frequency", fontsize=16)

# x-ticks major, -3, -1.5, 0, 1.5 and 3 only
axes.set_xticks([-3, -1.5, 0, 1.5, 3])

axes.set_title("")
# no spines top right and top left
axes.spines["top"].set_visible(False)
axes.spines["right"].set_visible(False)

axes.legend(frameon=False, fontsize=20)

plt.savefig(
    "/work/mh0033/m300883/Tel_MMLE/docs/defense/NAO_full_signal.png",
    dpi=500,
    bbox_inches="tight",
)

# %%

# %%

# %%
# New plot: histogram of mean-removed (centered) data
fig_centered, ax_centered = plt.subplots(1, 1, figsize=(34 * cm, 19 * cm))

ppc_first_ano_centered = ppc_first_ano - ppc_first_ano.mean()
ppc_last_ano_centered = ppc_last_ano - ppc_last_ano.mean()

ppc_first_ano_centered.plot.hist(
    ax=ax_centered,
    alpha=0.7,
    bins=np.arange(-4, 4.1, 0.5),
    color="k",
    label="1850-1859",
    edgecolor="white",
    linewidth=0.5,
)
ppc_last_ano_centered.plot.hist(
    ax=ax_centered,
    alpha=0.7,
    bins=np.arange(-4, 4.1, 0.5),
    color="r",
    label="2090-2099",
    edgecolor="white",
    linewidth=0.5,
)

# vline at zero (mean-removed)
ax_centered.axvline(0, color="k", linestyle="dashed", linewidth=3)
# ax_centered.axvline(
#     ppc_last_ano_centered.mean(), color="r", linestyle="dashed", linewidth=3
# )

ax_centered.tick_params(axis="both", which="major", labelsize=14)
ax_centered.set_xlabel("standarized NAO index (mean removed)", fontsize=20)
ax_centered.set_ylabel("Frequency", fontsize=16)
ax_centered.set_xticks([-3, -1.5, 0, 1.5, 3])
ax_centered.set_title("")
ax_centered.spines["top"].set_visible(False)
ax_centered.spines["right"].set_visible(False)
ax_centered.legend(frameon=False, fontsize=20)

plt.savefig(
    "/work/mh0033/m300883/Tel_MMLE/docs/defense/NAO_full_signal_centered.png",
    dpi=500,
    bbox_inches="tight",
)
# %%

# New plot: histogram with hatches showing the difference where 2090-2099 > 1850-1859 (mean-removed data)
fig_hatch, ax_hatch = plt.subplots(1, 1, figsize=(34 * cm, 19 * cm))

bins = np.arange(-4, 4.1, 0.5)
# Use centered (mean-removed) data
ppc_first_ano_centered = ppc_first_ano - ppc_first_ano.mean()
ppc_last_ano_centered = ppc_last_ano - ppc_last_ano.mean()
hist_first, _ = np.histogram(ppc_first_ano_centered, bins=bins)
hist_last, _ = np.histogram(ppc_last_ano_centered, bins=bins)
bin_centers = (bins[:-1] + bins[1:]) / 2
width = bins[1] - bins[0]

# Plot black bars (1850-1859)
bars_first = ax_hatch.bar(
    bin_centers,
    hist_first,
    width=width,
    alpha=0.7,
    color="k",
    label="1850-1859",
    edgecolor="white",
    linewidth=0.5,
    zorder=2,
)

# Plot red bars (2090-2099)
bars_last = ax_hatch.bar(
    bin_centers,
    hist_last,
    width=width,
    alpha=0.7,
    color="r",
    label="2090-2099",
    edgecolor="white",
    linewidth=0.5,
    zorder=3,
)

# Add hatched difference where 2090-2099 > 1850-1859
increase_label_added = False
for i, (center, val_first, val_last) in enumerate(
    zip(bin_centers, hist_first, hist_last)
):
    if val_last > val_first and (center > 1.5 or center < -1.5):
        ax_hatch.bar(
            center,
            val_last - val_first,
            width=width,
            bottom=val_first,
            color="none",
            hatch="//",
            linewidth=1.2,
            zorder=4,
            label="Increase" if not increase_label_added else None,
        )
        increase_label_added = True

# vline at zero (mean-removed)
# ax_hatch.axvline(
#     ppc_last_ano_centered.mean(), color="r", linestyle="dashed", linewidth=3
# )

ax_hatch.tick_params(axis="both", which="major", labelsize=14)
ax_hatch.set_xlabel("standarized NAO index (mean removed)", fontsize=20)
ax_hatch.set_ylabel("Frequency", fontsize=16)
ax_hatch.set_xticks([-3, -1.5, 0, 1.5, 3])
ax_hatch.set_title("")
ax_hatch.spines["top"].set_visible(False)
ax_hatch.spines["right"].set_visible(False)
ax_hatch.legend(frameon=False, fontsize=20)
ax_hatch.axvline(0, color="k", linestyle="dashed", linewidth=3, zorder=100)

plt.savefig(
    "/work/mh0033/m300883/Tel_MMLE/docs/defense/NAO_full_signal_hatched_centered.png",
    dpi=500,
    bbox_inches="tight",
)
# %%
