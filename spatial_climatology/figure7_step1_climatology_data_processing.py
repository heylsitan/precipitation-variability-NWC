import xarray as xr
import pandas as pd
import numpy as np

nw_file_path = '/data/nw_prev/mean_prev_1982_2020.nc'
ne_file_path = '/data/ne_prev/mean_prev_1982_2020.nc'
sw_file_path = '/data/sw_prev/mean_prev_1982_2020.nc'
se_file_path = '/data/se_prev/mean_prev_1982_2020.nc'
output_file = ''

ds_nw = xr.open_dataset(nw_file_path)
ds_ne = xr.open_dataset(ne_file_path)
ds_sw = xr.open_dataset(sw_file_path)
ds_se = xr.open_dataset(se_file_path)
nw_prev = ds_nw['tp']
ne_prev = ds_ne['tp']
sw_prev = ds_sw['tp']
se_prev = ds_se['tp']
total_prev = nw_prev + ne_prev + sw_prev + se_prev

con = xr.where(total_prev == 0, np.nan, sw_prev / total_prev)

contributition = xr.Dataset({'con':con})
contributition.to_netcdf(output_file)