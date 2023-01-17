import proplot as pplt
import matplotlib.pyplot as plt


def plot_scatter(ext_counts, tsurf, axes):
    extr_types = ["pos", "neg"]  # rows
    modes = ["NAO", "EA"]  # cols

    for i, extr_type in enumerate(extr_types):
        for j, mode in enumerate(modes):
            axes[i, j].scatter(
                x=tsurf,
                y=ext_counts.sel(extr_type=extr_type, mode=mode),
                alpha=0.5,
            )
    return axes


def extCount_tsurf_scatter(extc_tsuf_pairs, labels):
    """
    rows: pos/neg
    cols: NAO/EA
    scatter: extreme_count v.s surface temperature
    hue: different  dataset
    """
    fig, axes = pplt.subplots(nrows=2, ncols=2, figwidth=8, span=False, share=False)

    axes.format(
        abc="a",
        abcloc="ul",
        xlim=(-1, 6),
        suptitle=f"extreme counts v.s surface temperature in decadal time scale",
        xlabel="temperature (K)",
        ylabel="extreme count",
        grid=False,
        toplabels=["NAO", "EA"],
        leftlabels=["pos", "neg"],
    )

    for i, (ext_counts, t_surf) in enumerate(extc_tsuf_pairs):
        hs = plot_scatter(ext_counts, t_surf, label=labels[i], axes=axes)
    axes[1,1].legend(hs, ncols=1)
