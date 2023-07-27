#%%
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import proplot as pplt

# %%
def read_extrc(month, model="MPI_GE_onepct"):
    dir = f"/work/mh0033/m300883/Tel_MMLE/data/{model}/extreme_count/"
    filename = f"plev_50000_decade_first_{month}_extre_counts.nc"
    ds = xr.open_dataset(dir + filename).sel(confidence="true")
    ds["time"] = ds["time.year"]
    return ds.pc


def read_extrc_month(model="MPI_GE_onepct"):
    months = [
        "Jan",
        "Feb",
        "Mar",
        "Apr",
        "May",
        "Jun",
        "Jul",
        "Aug",
        "Sep",
        "Oct",
        "Nov",
        "Dec",
    ]
    mon_idx = xr.IndexVariable("month", months)
    ds = [read_extrc(month, model) for month in months]
    ds = xr.concat(ds, dim=mon_idx)
    return ds


#%%
# MPI_GE_onepct
MPI_GE_onepct = read_extrc_month("MPI_GE_onepct")
MPI_GE_onepct = MPI_GE_onepct.isel(time=[0, -1])
# %%
# MMLEA
models = ["MPI_GE", "CanESM2", "CESM1_CAM5", "MK36", "GFDL_CM3"]
models_ind = xr.IndexVariable("model", models)
months = [
    "Jan",
    "Feb",
    "Mar",
    "Apr",
    "May",
    "Jun",
    "Jul",
    "Aug",
    "Sep",
    "Oct",
    "Nov",
    "Dec",
]

MMLEA = [read_extrc_month(model) for model in models]
MMLEA = xr.concat(MMLEA, dim=models_ind)
MMLEA = MMLEA.sel(time=slice("1950-01-01", "2091-01-01"))

#%%
# Calculate the difference between the first and last decade
MMLEA_diff = MMLEA.isel(time=-1) - MMLEA.isel(time=0)

#%%
# Calculate the slope of the linear regression along the 'time' dimension
MMLEA_slope = MMLEA.polyfit(dim="time", deg=1).polyfit_coefficients.sel(degree=1)

# %%
def bar_hatch():
    pass

gs = pplt.GridSpec(nrows=1, ncols=2)
fig = pplt.figure(refwidth=2.2, span=False, share="labels")
# the right order of the models
models = ["MPI_GE", "CanESM2", "CESM1_CAM5", "MK36", "GFDL_CM3"]
colors_model = ["C1", "tab:blue", "tab:purple", "tab:cyan", "C4"]
model_color = dict(zip(models, colors_model))
months = [
    "Nov",
    "Oct",
    "Sep",
    "Aug",
    "Jul",
    "Jun",
    "May",
    "Apr",
    "Mar",
    "Feb",
    "Jan",
    "Dec",
]


for r, mode in enumerate(["NAO"]):
    for c, extr_type in enumerate(["pos", "neg"]):
        ax = fig.subplot(gs[r, c])

        for i, month in enumerate(months):
            data_month = MMLEA_slope.sel(mode=mode, extr_type=extr_type, month=month)
            data_cum = data_month.cumsum(axis=0).values
            for j, model in enumerate(models):
                data_model = data_month.sel(model=model).values
                bar_width = data_model
                starts = data_cum[j] - bar_width
                hatch = None
                zorder = 0
                if j > 0 and abs(data_cum[j]) < max(abs(data_cum[:j])):
                    hatch = "/////"
                    zorder = 10
                bar = ax.barh(
                    i,
                    bar_width,
                    fc=colors_model[j],
                    edgecolor="blue9",
                    label=model,
                    left=starts,
                    hatch = hatch,
                    width = 1.6,
                    alpha = 1,
                    zorder = zorder,
                )

                ax.format(
                    title=f"{mode} {extr_type}",
                    xminorticks="null",
                    yminorlocator="null",
                    grid=False,
                    xlim=(-0.4, 0.6),
                    ylim=(12, -1),
                    yticks=np.arange(11, -1, -1),
                    yticklabels=months,
                )

fig.legend(bar, ncols=5, loc="b", frame=False)


# fig.savefig("/work/mh0033/m300883/Tel_MMLE/docs/source/plots/monthly/month_extrc_diagram.png")


# %%
def bar_diff():

    gs = pplt.GridSpec(nrows=1, ncols=2)
    fig = pplt.figure(refwidth=2.2, span=False, share="labels")
    # the right order of the models
    models = ["MPI_GE", "CanESM2", "CESM1_CAM5", "MK36", "GFDL_CM3"]
    colors_model = ["C1", "tab:blue", "tab:purple", "tab:cyan", "C4"]
    model_color = dict(zip(models, colors_model))
    months = [
        "Nov",
        "Oct",
        "Sep",
        "Aug",
        "Jul",
        "Jun",
        "May",
        "Apr",
        "Mar",
        "Feb",
        "Jan",
        "Dec",
    ]



    for r, mode in enumerate(["NAO"]):
        for c, extr_type in enumerate(["pos", "neg"]):
            ax = fig.subplot(gs[r, c])

            for i, month in enumerate(months):
                data_month = MMLEA_slope.sel(mode=mode, extr_type=extr_type, month=month)
                # sort data_month by its value
                data_month = data_month.sortby(data_month)
                data_cum = data_month.cumsum(axis=0).values
                for j, model in enumerate(data_month.model.values):
                    data_model = data_month.sel(model=model).values
                    bar_width = data_model
                    start = data_cum[j] - bar_width
                    hatch = None
                    if j > 0 and abs(data_cum[j]) < abs(data_cum[j - 1]):
                        hatch = "/////"
                    bar = ax.barh(
                        i,
                        bar_width,
                        fc=model_color[model],
                        edgecolor="blue9",
                        label=model,
                        left=start,
                        hatch = hatch,
                        width = 1.6,

                    )

                    ax.format(
                        title=f"{mode} {extr_type}",
                        xminorticks="null",
                        yminorlocator="null",
                        grid=False,
                        xlim=(-0.4, 0.6),
                        ylim=(12, -1),
                        yticks=np.arange(11, -1, -1),
                        yticklabels=months,
                    )

    fig.legend(bar, ncols=5, loc="b", frame=False)

# %%
