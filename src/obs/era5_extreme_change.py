# %%
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from PyEMD import EEMD


# %%
def nao_without_decade(nao, method="EEMD", combine_num=2, rolling_years=[10, 20, 30]):
    NEW_nao = []
    for month in [6, 7, 8]:
        month_nao = nao.sel(time=nao.time.dt.month == month)

        if method == "EEMD":
            new_nao, eIMFs = remove_decade_eemd(month_nao, combine_num)
        elif method == "running_mean":
            new_nao = remove_decade_running_mean(month_nao, rolling_years)
        NEW_nao.append(new_nao)

    NEW_nao = xr.concat(NEW_nao, dim="time")  # concat the 3 months
    NEW_nao = NEW_nao.sortby("time")
    NEW_nao = NEW_nao.expand_dims({"mode": ["NAO"]})
    return NEW_nao


def remove_decade_running_mean(month_nao, rolling_years=[10, 20, 30]):
    rolling_mean = month_nao
    Rollings = []
    for year in rolling_years:
        rolling_mean = rolling_mean.rolling(
            time=year, center=True
        ).mean()  # do rolling on the rolling_mean before.
        rolling_mean = rolling_mean.interpolate_na(
            dim="time", method="slinear", fill_value="extrapolate"
        )
        Rollings.append(rolling_mean)

    rolling_mean = xr.concat(Rollings, dim="rolling_years")
    rolling_mean = rolling_mean.sum(dim="rolling_years")
    new_month_nao = month_nao - rolling_mean
    return new_month_nao


def remove_decade_eemd(month_nao, combine_num=2):
    # fix the random seed for reproducibility
    np.random.seed(42)
    eemd = EEMD()
    eIMFs = eemd.eemd(month_nao.values)
    new_nao = eIMFs[:combine_num].sum(axis=0)  # the first 3 IMFs
    new_nao = month_nao.copy(data=new_nao)
    return new_nao, eIMFs


def count_extreme(nao):
    nao_first = nao.sel(time=slice("1941", "1981")).values.flat
    nao_last = nao.sel(time=slice("1982", "2022")).values.flat
    first_pos = np.sum(nao_first > 1.5)
    first_neg = np.sum(nao_first < -1.5)
    last_pos = np.sum(nao_last > 1.5)
    last_neg = np.sum(nao_last < -1.5)
    return first_pos, first_neg, last_pos, last_neg


def plot_era_nao_index(nao, NEW_nao, ax):
    nao = nao.sel(time=slice("1941", "2022"))
    NEW_nao = NEW_nao.sel(time=slice("1941", "2022"))

    first_pos_org, first_neg_org, last_pos_org, last_neg_org = count_extreme(nao)
    first_pos_new, first_neg_new, last_pos_new, last_neg_new = count_extreme(NEW_nao)

    nao.plot.line(x="time", color="gray", alpha=0.8, ax=ax, lw=1)
    NEW_nao.plot.line(x="time", color="black", alpha=0.8, ax=ax, lw=1)

    # vline at x = 1981
    xmin, xmax = ax.get_xlim()
    xmid = (xmin + xmax) / 2
    ax.axvline(x=xmid, color="g", linestyle="--")
    ax.set_yticks([-3, -1.5, 0, 1.5, 3])

    # hline at y= 1.5 and -1.5
    ax.axhline(y=1.5, color="g", linestyle="--")
    ax.axhline(y=-1.5, color="g", linestyle="--")

    # put the count as text on the plot
    # pos
    ax.text(
        0.20,
        0.99,
        f"1941-1981: {first_pos_org}",
        transform=ax.transAxes,
        fontsize=7,
        color="gray",
        alpha=0.8,
        verticalalignment="top",
    )
    ax.text(
        0.55,
        0.99,
        f"1982-2022: {last_pos_org}",
        transform=ax.transAxes,
        fontsize=7,
        color="gray",
        alpha=0.8,
        verticalalignment="top",
    )

    ax.text(
        0.20,
        0.93,
        f"1941-1981: {first_pos_new}",
        transform=ax.transAxes,
        fontsize=7,
        color="black",
        alpha=0.8,
        verticalalignment="top",
    )
    ax.text(
        0.55,
        0.93,
        f"1982-2022: {last_pos_new}",
        transform=ax.transAxes,
        fontsize=7,
        color="black",
        alpha=0.8,
        verticalalignment="top",
    )

    # neg
    ax.text(
        0.20,
        0.15,
        f"1941-1981: {first_neg_org}",
        transform=ax.transAxes,
        fontsize=7,
        color="gray",
        alpha=0.8,
        verticalalignment="top",
    )
    ax.text(
        0.55,
        0.15,
        f"1982-2022: {last_neg_org}",
        transform=ax.transAxes,
        fontsize=7,
        color="gray",
        alpha=0.8,
        verticalalignment="top",
    )

    ax.text(
        0.20,
        0.10,
        f"1941-1981: {first_neg_new}",
        transform=ax.transAxes,
        fontsize=7,
        color="black",
        alpha=0.8,
        verticalalignment="top",
    )
    ax.text(
        0.55,
        0.10,
        f"1982-2022: {last_neg_new}",
        transform=ax.transAxes,
        fontsize=7,
        color="black",
        alpha=0.8,
        verticalalignment="top",
    )

    return ax
