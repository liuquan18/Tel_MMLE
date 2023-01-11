# Analysis process

## Depedent Index from dynamical spatial patterns
The indexes here are used to show the main plots in the paper.
- **Vertically**, *Horizontal-vertical EOF* (dependent) takes the whole troposphere as a whole, the magnitude of the spatial patterns strengthen, the shift of the centers of actions are mild (compared to horizontal EOF).
- **Temporally**, *dynamical spatial patterns* makes sure that the index also reveals the change of the spatial patterns. 

Since the dynamical spatial patterns are adopted, two ways of standardization are applied here 

![three series](plots/first10_last10/standard_ways.png)

### 1. standardize the index of first10 years and last10 years seperately
note that here the index of first10 years and last10 years are generated from different spaital patterns. `first10_first` and `last10_last` are standarded with the mean and std of its own (specificlly the mean and std of `all_first` and `all_last`.)

- generate `all_first` and `all_last` and standard them with their own means and stds. 
    
    ```python
    # standard levels       ---- make  all levels standard.
    # dependent             ---- the chosen vetical strategy.
    # fixed pattern         ---- fix the pattern to be `first` or `last`.
    # standard pc           ---- standard the pc with its own mean and std
    _,all_first,_ = season_eof(geopotential_height,nmode=2,
    window=10,fixed_pattern='first',independent = False,standard=True)

    _,all_last,_ = season_eof(geopotential_height,nmode=2,
    window=10,fixed_pattern='first',independent = False,standard=True)
    ```

- select the first10 (last10) years from all_first name as `first10_first ` (`last10_last`)

    ```python
    first10_first = all_first.isel(time = slice(0,10))

    last10_last = all_last.isel(time = slice(-10,len(all_last)))
    ```


### 2. Standardize the index from different patterns with mean and std of index from one common pattern (all-patterns).
 `all_first`, `all_last`, and `all_all` are generated. then `first10_first` and `last10_last` are selected. both of them are standard with the temporal mean and std of `all_all`.

- generate the index `all_first`,`all_last` and `all_all`, with the pc *unstandardized*.
    ```python
    _,all_first,_ = season_eof(geopotential_height,nmode=2,
    window=10,fixed_pattern='first',independent = False,standard=False)

    _,all_last,_ = season_eof(geopotential_height,nmode=2,
    window=10,fixed_pattern='first',independent = False,standard=False)

    _,all_all,_ = season_eof(geopotential_height,nmode=2,
    window=10,fixed_pattern='all',independent = False,standard=False)    
    ```

- select the first10 years and last10 years
    ```python
    first10_first = all_first.sel(time = slice(0,10))

    last10_last = all_last.sel(time = slice(-10,len(all_last)))
    ```

- standard them with the mean and std of `all_all`.
    ```python
    first10_first = Normalize(first10_first,all_all)
    last10_last = Normalize(last10_last, all_all)
    ```


## PDF of index from dynamic spatial patterns
just use the xarray function
```python
    xr.DataArray.plot.hist()
```

## Extreme count of index
get the vertical profile of the extreme count
```python
first10_ext_count = period_extreme_count(first10_first,dim = ("time","ens"),threshold = 2, standard = False)

last10_ext_count= period_extreme_count(first10_last,dim = ("time","ens"),threshold = 2, standard = False)
```
