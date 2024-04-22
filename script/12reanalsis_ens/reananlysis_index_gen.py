# %%
import numpy as np

import src.MMLE_TEL.index_generator as index_generate


import proplot as pplt
import src.plots.statistical_overview as stat_overview
import matplotlib.pyplot as plt


# %%
from src.reanalysis.decompose import EOF_reanalysis


# %%
# ERA5, 1950 - 2022, 40 years
ERA5_40 = EOF_reanalysis("ERA5", group_size=40, start_year="1950", end_year="2022")
ERA5_40.save_eof()
ERA5_40.plot_spatial_index(save=True)

# %%
# RA5, 1979 - 2022, 20 years
ERA5_20 = EOF_reanalysis("ERA5", group_size=20, start_year="1979", end_year="2022")
ERA5_20.save_eof()
ERA5_20.plot_spatial_index(save=True)

# %%
# RA5_allens, 1950 - 2022, 40 years
ERA5_allens_40 = EOF_reanalysis(
    "ERA5_allens", group_size=40, start_year="1950", end_year="2022"
)
ERA5_allens_40.save_eof()
ERA5_allens_40.plot_spatial_index(save=True)
# %%
# CR20, 1950 - 2015, 40 years
CR20_40 = EOF_reanalysis("CR20", group_size=40, start_year="1950", end_year="2015")
CR20_40.save_eof()
CR20_40.plot_spatial_index(save=True)

# %%
# CR20_allens, 1950 - 2015, 40 years
CR20_allens_40 = EOF_reanalysis(
    "CR20_allens",
    group_size=40,
    start_year="1850",
    end_year="2015",
    external_forcing="quadratic_trend",
    period_dec=False,
)
CR20_allens_40.save_eof()
CR20_allens_40.plot_spatial_index(save=True)
# %%
# period decomposition
CR20_allens_40 = EOF_reanalysis(
    "CR20_allens",
    group_size=40,
    start_year="1850",
    end_year="2015",
    external_forcing="quadratic_trend",
    period_dec=True,
)

# %%
CR20_allens_40.first_eof_std.to_netcdf(
    f"/work/mh0033/m300883/Tel_MMLE/data/CR20_allens/EOF_result/periodic_first_40_eof_std.nc"
)
CR20_allens_40.last_eof_std.to_netcdf(
    f"/work/mh0033/m300883/Tel_MMLE/data/CR20_allens/EOF_result/periodic_last_40_eof_std.nc"
)

# %%
