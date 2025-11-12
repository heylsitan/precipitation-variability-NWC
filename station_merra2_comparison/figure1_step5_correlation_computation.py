import numpy as np
import xarray as xr

years = np.arange(1982,2021)
for year in years:
    print(year)
    input_file = f""
    ds = xr.open_dataset(input_file)
    ds_mean = ds.std(dim='time', skipna=False)  #skipna=True is used to skip nan value calculations
    output_file = f""
    ds_mean.to_netcdf(output_file)