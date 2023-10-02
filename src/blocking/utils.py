#%%
import xarray as xr
import numpy as np

#%%

def reg_lens(state_arr):
    sa=np.asarray(state_arr)
    n=len(sa)
    
    if n==0:
        return (None,None)
    
    else:
        #we create a truth array for if a state is about to change:
        x=np.array(sa[1:]!=sa[:-1])
        #we create an array containing those entries of sa where x is true
        #and append the final entry:
        y=np.append(np.where(x),n-1)
        #we create an array of persistent state lengths:
        L=np.diff(np.append(-1,y))
        #and an array of those state values:
        S=sa[y]
        return (L,S)


#%%
def persistent_state(ix,state,min_pers):
    L,S=reg_lens(ix)

    keepS=S==state
    keepL=L>min_pers
    keep_points=(np.repeat(keepL,L)*np.repeat(keepS,L))

    return keep_points