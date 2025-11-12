import xarray as xr
import numpy as np
import pandas as pd

# Used to calculate precipitation for each source region
date_list = pd.date_range('19820101', '20201231')
summer_dates = [
    date for date in date_list
    if date.month in [6, 7, 8] and not (date.month == 6 and date.day == 1)
]

for date in summer_dates:
    print(date)
    input_file = f"" 
    output_file = f''
    ds = xr.open_dataset(input_file)
    tp = ds['p_track_lower'] + ds['p_track_upper']
    tp = xr.Dataset({'tp':tp})
    tp.to_netcdf(output_file)