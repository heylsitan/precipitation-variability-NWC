import xarray as xr
import numpy as np
import pandas as pd

years = np.arange(1982,2021)
for year in years:
    print(year)
    prr_file = f''
    tp_file = f''
    ds_prr = xr.open_dataset(prr_file)
    ds_tp = xr.open_dataset(tp_file)
    prr_value = ds_prr['prr']
    tp_value = ds_tp['tp']
    intrenal_tp = tp_value*prr_value
    in_tp = xr.Dataset({'tp':intrenal_tp})
    output_in_file = f''
    in_tp.to_netcdf(output_in_file)