import xarray as xr
import numpy as np
import os

# Folder path
folder_path = '/data/merra2_prev_nwc'

dry_years = [1983, 1985, 1986, 1994, 1995]
# Get paths for all .nc files
dry_file_paths = [os.path.join(folder_path, f'prev_{year}.nc') for year in dry_years]

# Initialize a list to store data
dry_data_list = []

# Read each .nc file and add to list
for file_path in dry_file_paths:
    ds = xr.open_dataset(file_path)  # Open .nc file
    dry_data_list.append(ds['tp'])  # Assume variable name is 'precipitation', adjust according to actual situation

# Combine all data into one dataset
dry_combined_data = xr.concat(dry_data_list, dim='time')

# Calculate average of all data
dry_mean_data = dry_combined_data.mean(dim='time')


wet_years = [1983, 1985, 1986, 1994, 1995]
# Get paths for all .nc files
wet_file_paths = [os.path.join(folder_path, f'prev_{year}.nc') for year in wet_years]

# Initialize a list to store data
wet_data_list = []

# Read each .nc file and add to list
for file_path in wet_file_paths:
    ds = xr.open_dataset(file_path)  # Open .nc file
    wet_data_list.append(ds['tp'])  # Assume variable name is 'precipitation', adjust according to actual situation

# Combine all data into one dataset
wet_combined_data = xr.concat(wet_data_list, dim='time')

# Calculate average of all data
wet_mean_data = wet_combined_data.mean(dim='time')

difference_mean = wet_mean_data - dry_mean_data

difference_mean.to_netcdf()