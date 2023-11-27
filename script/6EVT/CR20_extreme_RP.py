#%%
from pyextremes import get_extremes
from pyextremes.plotting import plot_extremes
from pyextremes import get_return_periods
from pyextremes import EVA

import xarray as xr
import numpy as np
import pandas as pd
# %%
EOF = xr.open_dataset("/work/mh0033/m300883/Tel_MMLE/data/CR20/EOF_result/all_40_eof.nc")
PC = EOF.pc
# %%
NAO = PC.sel(mode = 'NAO')
# %%
NAO = NAO.squeeze()
df_NAO = NAO.to_dataframe()
# %%
data = df_NAO['pc']
#%%
model = EVA(data)

#%%
extremes = model.get_extremes("BM", extremes_type="low", block_size = '10Y')
#%%
model.plot_extremes()

#%%
model.fit_model()
#%%
model.model.get_return_value(model.model.pdf(1.5))



#%%
#%%
summary = model.get_summary(
    return_period=[1, 2, 5, 10, 25, 50, 100, 250, 500, 1000],
    alpha=0.95,
    n_samples=1000,
)
#%%
print(summary)
model.plot_diagnostic(alpha=0.95)

# %%
pos_extremes = get_extremes(data, "BM", extremes_type="high",block_size = '10Y')
# %%
neg_extremes = get_extremes(data, "BM", extremes_type="low",block_size = '10Y')


# %%
pos_return_periods = get_return_periods(
    ts=data,
    extremes=pos_extremes,
    extremes_method="BM",
    extremes_type="high",
    block_size="10Y",
    return_period_size="365.2425D",
    plotting_position="weibull",
)
pos_return_periods.sort_values("return period", ascending=False).head()
# %%
