#%%
from mdutils.mdutils import MdUtils

#%%
def doc_quick_plots(file_name, title, plot_dir, prefix):
    qp = MdUtils(file_name, title="")  # qick plots,no title here
    qp.new_header(level = 1, title = title)

    # overview
    qp.new_header(level=2, title="statistical overview")
    qp.new_paragraph("The spatial pattern and distribution of NAO and EA")
    qp.new_line(
        qp.new_inline_image(
            text="spatial_patternv",
            path=plot_dir + prefix + "spatial_pattern_violin500hpa.png",
        )
    )
    qp.new_paragraph("The spatial pattern and distribution of NAO and EA")
    qp.new_line(
        qp.new_inline_image(
            text="spatial_patternh",
            path=plot_dir + prefix + "spatial_pattern_hist500hpa.png",
        )
    )

    qp.new_paragraph("the violin plot of NAO and EA for all vertical levels")
    qp.new_line(
        qp.new_inline_image(
            text="violin profile", path=plot_dir + prefix + "violin_profile.png"
        )
    )

    # extreme count
    qp.new_header(level=2, title="extreme count")
    qp.new_paragraph("the extreme count of NAO and EA for all vertical levels")

    qp.new_header(level=3, title="NAO profile")
    qp.new_line(
        qp.new_inline_image(
            text="NAO_extreme_count_profile",
            path=plot_dir + prefix + "NAO_extreme_count_profile.png",
        )
    )

    qp.new_header(level=3, title="EA profile")
    qp.new_line(
        qp.new_inline_image(
            text="EA_extreme count profile",
            path=plot_dir + prefix + "EA_extreme_count_profile.png",
        )
    )

    # return period
    qp.new_header(level=2, title="return period")

    qp.new_paragraph("the return period of NAO and EA at 500hpa in different periods")
    qp.new_header(level=3, title="500hpa scatter")
    qp.new_line(
        qp.new_inline_image(
            text="500hpa scatter return period",
            path=plot_dir + prefix + "NAO_return_period_scatter.png",
        )
    )

    qp.new_header(level=3, title="NAO profile")
    qp.new_line(
        qp.new_inline_image(
            text="NAO profile return period",
            path=plot_dir + prefix + "NAO_return_period_profile.png",
        )
    )

    qp.new_header(level=3, title="EA profile")
    qp.new_line(
        qp.new_inline_image(
            text="EA profile return period",
            path=plot_dir + prefix + "EA_return_period_profile.png",
        )
    )

    qp.new_header(level=2, title="extreme spatial patterns")
    qp.new_line(
        qp.new_inline_image(
            text="extreme sptial pattern",
            path=plot_dir + prefix + "extreme_spatial_pattern_1000hpa.png",
        )
    )

    qp.new_header(level=2, title="Influence on surface temperature")

    qp.new_header(level=3, title="NAO index at 500hpa")
    qp.new_line(
        qp.new_inline_image(
            text="tsurf", path=plot_dir + prefix + "composite_tsurf_NAO.png"
        )
    )

    qp.new_header(level=3, title="EA index at 500hpa")
    qp.new_line(
        qp.new_inline_image(
            text="tsurf", path=plot_dir + prefix + "composite_tsurf_EA.png"
        )
    )

    # Create a table of contents
    # qp.new_table_of_contents(table_title="Contents", depth=2)
    qp.create_md_file()


# %%
