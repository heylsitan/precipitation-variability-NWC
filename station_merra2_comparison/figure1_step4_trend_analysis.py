import numpy as np
import xarray as xr
import pandas as pd
import geopandas as geo
import salem
from pathlib import Path

years = np.arange(1980,2021)
input_folder = ''
output_folder = ''
china_shp = geo.read_file('/home/northwest/northwest_China.shp')

for year in years:
    print(year)
    file_name = f"total_precipitation_{year}.nc"  # Corrected line
    input_path = Path(input_folder) / file_name
    output_file = Path(output_folder) / file_name
    total_precipitation = xr.open_dataset(input_path)
    tp = total_precipitation['tp']
    china_tp = tp.salem.roi(shape=china_shp)
    china_tp.to_netcdf(output_file)