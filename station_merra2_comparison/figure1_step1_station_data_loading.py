import glob
import numpy as np
import xarray as xr
import pandas as pd

years = np.arange(1982,2021)
datelist = pd.date_range('19820101','20201231')
summer_dates = [date for date in datelist if date.month in [6, 7, 8]]

for date in summer_dates:
    print(date.strftime('%Y%m%d'))
    # ====== Read data ======
    pattern = f""
    files = glob.glob(pattern)
    ds = xr.open_dataset(files[0])
    merra2_tp = ds['PRECTOT']*3.6
    merra2_tp["time"] = merra2_tp["time"] - np.timedelta64(30, "m")
    
    tp = xr.Dataset({'tp': merra2_tp})
    tp.to_netcdf(f"")