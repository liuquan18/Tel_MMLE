import proplot as pplt
import src.EVT.return_period as EVT


def return_period_scatter(pos,mpos,neg,mneg, mode, periods, labels,hlayers=50000):
    """
    return period plot
    """
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

    all_type = [pos,neg]
    media_type = [mpos,mneg]
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c']
    markers = ["+","*","d"]

    for i, extr_type in enumerate(all_type):
        for j, period in enumerate(extr_type):
            # all points
            axes[i].scatter(
                x = "return period",
                y = "pc",
                data = period,
                label = labels[j],
                color = colors[j],
            )

            # media points
            axes[i].scatter(
                x = "return period",
                y = "pc",
                data = media_type[i][j],
                color = "k",
                marker = markers[j],
                label = labels[j]+'med',
                s = 80,
                zorder = 100+j
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
        ix.scatter(x="return period", y="pc", data=neg[-1], color=colors[-1])
        ix.set_xlabel(None)
        ix.set_ylabel(None)

    axes[0].legend(loc="lr", ncols=1)
    axes[0].format(title="pos")
    axes[1].format(title="neg")


def return_period_profile(pos, neg, index, mode,labels):
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
    extrs = ['pos','neg']
    extr_data = [pos,neg]
    styles = ['-','--','dotted']
    for i, extr_type in enumerate(extrs):
        for j in range(pos.shape[-1]):
            axes[i].plot(extr_data[i][:,j],y,color = 'k',linestyle = styles[j],label = labels[j])
            axes[i].set_title(extr_type)

    axes[1].set_xlim(0, 5)
    axes[0].set_xlim(0, 5)

    axes[0].legend(loc="lr", ncols=1)
    axes[1].legend(loc="lr", ncols=1)

    for ax in axes:
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)

    axes[0].format(title="pos")
    axes[1].format(title="neg")
