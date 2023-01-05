import proplot as pplt
import src.EVT.return_period as EVT


def return_period_scatter(index, mode, hlayers=50000):
    """
    return period plot
    """

    (
        first10_all_pos,
        last10_all_pos,
        first10_median_pos,
        last10_median_pos,
        first10_all_neg,
        last10_all_neg,
        first10_median_neg,
        last10_median_neg,
    ) = EVT.mode_return_period(index, mode=mode, hlayers=hlayers)

    fig, axes = pplt.subplots(nrows=1, ncols=2, figwidth=8, span=False, share=False)

    axes.format(
        abc="a",
        abcloc="ul",
        xlim=(0, 10),
        xminorlocator="null",
        yminorlocator="null",
        suptitle=f"return period of {mode} index at {hlayers/100:.0f}hpa",
        xlabel="return period / yr",
        ylabel="pc / std",
        grid=False,
    )

    # pos
    axes[0].scatter(
        x="return period",
        y="pc",
        data=first10_all_pos,
        label="first10",
    )
    axes[0].scatter(
        x="return period",
        y="pc",
        data=last10_all_pos,
        color="r",
        label="last10",
    )

    # pos media
    axes[0].scatter(
        x="return period",
        y="pc",
        data=first10_median_pos,
        color="k",
        label="first10 median",
        marker="+",
        s = 80,
    )
    axes[0].scatter(
        x="return period",
        y="pc",
        data=last10_median_pos,
        color="k",
        label="last10 median",
        marker="*",
        s = 80,
    )

    # neg
    axes[1].scatter(
        x="return period",
        y="pc",
        data=first10_all_neg,
        label="first10",
    )
    axes[1].scatter(
        x="return period",
        y="pc",
        data=last10_all_neg,
        color="r",
        label="last10",
    )
    # neg media
    axes[1].scatter(
        x="return period",
        y="pc",
        data=first10_median_neg,
        color="k",
        label="first10 median",
        marker="+",
        s= 80
    )
    axes[1].scatter(
        x="return period",
        y="pc",
        data=last10_median_neg,
        color="k",
        label="last10 median",
        marker="*",
        s = 80
    )

    # a zoom out for NAO very extreme events.
    if mode == "NAO":
        ix = axes[1].inset_axes([6, -2.8, 3.8, 0.2], transform="data", zoom=False)
        ix.format(
            xlim=(0, 300),
            ylim=(-4, -2.6),
            xminorlocator="null",
            yminorlocator="null",
            title="zoom out",
        )
        ix.scatter(x="return period", y="pc", data=last10_all_neg, color="r")
        ix.set_xlabel(None)
        ix.set_ylabel(None)

    axes[0].legend(loc="lr", ncols=1)
    axes[0].format(title="pos")
    axes[1].format(title="neg")


def return_period_profile(pos, neg, index, mode):
    fig, axes = pplt.subplots(nrows=1, ncols=2, figwidth=8, span=False, share=False)

    axes.format(
        abc="a",
        abcloc="ul",
        xlim=(0, 5),
        xminorlocator="null",
        yminorlocator="null",
        suptitle=f"media return period of {mode} index",
        xlabel="return period / yr",
        ylabel="gph/hpa",
        ylim=(1000, 200),
        grid=False,
    )

    y = index.hlayers.values / 100
    axes[0].plot(pos[:, 0], y, "k-", label="first10 median")

    axes[0].plot(pos[:, 1], y, "k--", label="last10 median")

    axes[1].plot(neg[:, 0], y, "k-", label="first10 median")
    axes[1].plot(neg[:, 1], y, "k--", label="last10 median")

    axes[1].set_xlim(0, 5)
    axes[0].set_xlim(0, 5)

    axes[0].legend(loc="lr", ncols=1)
    axes[1].legend(loc="lr", ncols=1)

    for ax in axes:
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)

    axes[0].format(title="pos")
    axes[1].format(title="neg")
