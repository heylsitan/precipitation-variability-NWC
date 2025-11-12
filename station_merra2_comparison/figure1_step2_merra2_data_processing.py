import numpy as np
import pandas as pd 
import xarray as xr

date_list = pd.date_range("19820101", "20201231")
summer_dates = [date for date in date_list if date.month in [6, 7, 8]]

for summer_date in summer_dates:
    print(summer_date.strftime('%Y%m%d'))
    file_path = f''
    ds = xr.open_dataset(file_path)
    tp = ds['tp']
    first_time = tp['time'].isel(time=0)
    sum_tp = tp.sum(dim="time", skipna=False)
    ds_output = sum_tp.expand_dims(time=[first_time.values])
    ds_output.to_netcdf(f'')