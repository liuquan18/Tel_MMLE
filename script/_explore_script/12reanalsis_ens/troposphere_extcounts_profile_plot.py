#%%
import xarray as xr
import src.plots.extreme_plot as extreme_plot
import proplot as pplt
import numpy as np
import matplotlib.pyplot as plt
# %%
import importlib
importlib.reload(extreme_plot)


# %%
def extreme_count_rean(model, plevs):
    odir = f"/work/mh0033/m300883/Tel_MMLE/data/{model}/extreme_count/"

    First = []
    Last = []

    for plev in plevs:
        first = xr.open_dataset(f"{odir}first_plev{plev}_extc.nc")
        last = xr.open_dataset(f"{odir}last_plev{plev}_extc.nc")

        # a new dimension called plev
        first['plev'] = plev
        last['plev'] = plev

        First.append(first)
        Last.append(last)


    # concat along the plev dimension
    First = xr.concat(First, dim='plev')
    Last = xr.concat(Last, dim='plev')

    First = First.sortby('plev')
    Last = Last.sortby('plev')

    # divided by ensemble size
    if model == 'CR20_allens':
        first_count = First.pc/(4 * 80)
        last_count = Last.pc/(4 * 80)
    else:
        first_count = First.pc/(4)
        last_count = Last.pc/(4)
    return first_count, last_count

def extreme_count_MPI(model):
    vertical_eof = 'ind'
    prefix = f"/work/mh0033/m300883/Tel_MMLE/data/{model}/extreme_count/troposphere_{vertical_eof}_decade_first_JJA_"
    first = xr.open_dataset(f"{prefix}first_count.nc")
    last = xr.open_dataset(f"{prefix}last_count.nc")

    # divided by ensemble size
    first_count = first.pc/100
    last_count = last.pc/100
    return first_count, last_count


#%%
model = "CR20_allens"
plevs = [92500, 85000, 70000, 50000, 40000, 30000, 20000]
CR20_all_first, CR20_all_last = extreme_count_rean(model, plevs)


# %%
model = 'CR20'
plevs = [1000, 975, 950, 925, 900, 850, 800, 750, 700, 650, 600, 550, 500, 450, 400, 350, 300, 250, 200]
CR20_mean_first, CR20_mean_last = extreme_count_rean(model, plevs)

#%%
model = 'MPI_GE_onepct'
onepct_first, onepct_last = extreme_count_MPI(model)

#%%
model = 'MPI_GE'
MPI_first, MPI_last = extreme_count_MPI(model)

# %%
fig, axes = pplt.subplots(
    space=0,
    width=150 / 25.4,
    height=150 / 25.4,
    hspace=5,
    wspace=2,
    nrows=2,
    ncols=2,
    sharex=False,
    sharey=False,
)
axes.format(
    abc=True,
    abcloc="ul",
    abcstyle="a",
    xlabel="Extreme event count",
    ylabel="Pressure level (hPa)",
    ylim=(1000, 200),
    xminorticks="null",
    yminorticks="null",
    grid=False,
    toplabels=("pos", "neg"),
    # leftlabels=("NAO", "EA"),
    # xlocator=20,
)

ens_labels = ["first10", "last10"]
for i, count in enumerate([onepct_first, onepct_last]):
    extreme_plot.plot_extreme_count(
        count.sel(mode = 'NAO', extr_type ="pos"),
        axes[0, 0],
        label=ens_labels[i],
    )
    extreme_plot.plot_extreme_count(
        count.sel(mode = 'NAO', extr_type ="neg"),
        axes[0, 1],
        label=ens_labels[i],
    )

mean_labels = ["first40", "last40"]
for i, count in enumerate([CR20_mean_first, CR20_mean_last]):
    extreme_plot.plot_extreme_count(
        count.sel(mode = 'NAO', extr_type ="pos"),
        axes[1, 0],
        label=mean_labels[i],
    )
    extreme_plot.plot_extreme_count(
        count.sel(mode = 'NAO', extr_type ="neg"),
        axes[1, 1],
        label=mean_labels[i],
    )
axes[0,0].set_xlim(1,4)
axes[0,1].set_xlim(1,4)

axes[1,0].set_xlim(0,7)
axes[1,1].set_xlim(0,7)

# axes[0,0].set_xticks(np.arange(1,4.1,0.5))
# axes[0,1].set_xticks(np.arange(1,4.1,0.5))

# no yticks ylabel for the second columns
for ax in axes[:, 1]:
    ax.format(ylabel="")
plt.setp(axes[0,1].get_yticklabels(), visible=False)
plt.setp(axes[1,1].get_yticklabels(), visible=False)

axes[0, 0].legend(loc="lr", ncols=1, frame=True)
axes[1, 0].legend(loc="lr", ncols=1, frame=True)
for ax in axes:
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)


plt.savefig("/work/mh0033/m300883/Tel_MMLE/docs/source/plots/paper_main/profile.pdf")
# %%


#### supplementary for MPI_GE and CR20_allens

# %%
fig, axes = pplt.subplots(
    space=0,
    width=150 / 25.4,
    height=150 / 25.4,
    hspace=5,
    wspace=2,
    nrows=2,
    ncols=2,
    sharex=False,
    sharey=False,
)
axes.format(
    abc=True,
    abcloc="ul",
    abcstyle="a",
    xlabel="Extreme event count",
    ylabel="Pressure level (hPa)",
    ylim=(1000, 200),
    xminorticks="null",
    yminorticks="null",
    grid=False,
    toplabels=("pos", "neg"),
    # leftlabels=("NAO", "EA"),
    # xlocator=20,
)

ens_labels = ["first10", "last10"]
for i, count in enumerate([MPI_first, MPI_last]):
    extreme_plot.plot_extreme_count(
        count.sel(mode = 'NAO', extr_type ="pos"),
        axes[0, 0],
        label=ens_labels[i],
    )
    extreme_plot.plot_extreme_count(
        count.sel(mode = 'NAO', extr_type ="neg"),
        axes[0, 1],
        label=ens_labels[i],
    )

mean_labels = ["first40", "last40"]
for i, count in enumerate([CR20_all_first, CR20_all_last]):
    extreme_plot.plot_extreme_count(
        count.sel(mode = 'NAO', extr_type ="pos"),
        axes[1, 0],
        label=mean_labels[i],
    )
    extreme_plot.plot_extreme_count(
        count.sel(mode = 'NAO', extr_type ="neg"),
        axes[1, 1],
        label=mean_labels[i],
    )
axes[0,0].set_xlim(1,4)
axes[0,1].set_xlim(1,4)

axes[1,0].set_xlim(1,4)
axes[1,1].set_xlim(1,4)

# axes[0,0].set_xticks(np.arange(1,4.1,0.5))
# axes[0,1].set_xticks(np.arange(1,4.1,0.5))

# no yticks ylabel for the second columns
for ax in axes[:, 1]:
    ax.format(ylabel="")
plt.setp(axes[0,1].get_yticklabels(), visible=False)
plt.setp(axes[1,1].get_yticklabels(), visible=False)

axes[0, 0].legend(loc="lr", ncols=1, frame=True)
axes[1, 0].legend(loc="lr", ncols=1, frame=True)

for ax in axes:
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)

plt.savefig("/work/mh0033/m300883/Tel_MMLE/docs/source/plots/paper_supplymentary/MPI_GE_CR20_all_profile.pdf")

# %%
