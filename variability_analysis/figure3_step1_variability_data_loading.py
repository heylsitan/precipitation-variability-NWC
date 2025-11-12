import xarray as xr
import pandas as pd
import numpy as np

# Create a date range containing multiple years
years = np.arange(1982,2021)

for year in years:
    print(year)
    tag_tp_file_path = f'data/wam2/result.nc'
    tp_file_path =f'data/merra2/tp.nc'
    output_file =f''
    ds_e_track = xr.open_dataset(tag_tp_file_path)
    ds_tp = xr.open_dataset(tp_file_path)
    e_track = ds_e_track['tag_tp']
    tp = ds_tp['tp']
    prr_value = xr.where(tp == 0, 0, e_track / tp)
    prr = xr.Dataset({'prr':prr_value})
    prr.to_netcdf(output_file)