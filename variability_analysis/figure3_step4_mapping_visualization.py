import xarray as xr
import numpy as np

years = np.arange(1982,2021)

for year in years:
    print(year)
    nw_file_path = f"/data/track_tp/nw_year/tp_{year}.nc"
    ne_file_path = f"/data/track_tp/ne_year/tp_{year}.nc"
    sw_file_path = f"/data/track_tp/sw_year/tp_{year}.nc"
    se_file_path = f"/data/track_tp/se_year/tp_{year}.nc"

    # Read precipitation for each year
    nw_tp = xr.open_dataset(nw_file_path)
    ne_tp = xr.open_dataset(ne_file_path)
    sw_tp = xr.open_dataset(sw_file_path)
    se_tp = xr.open_dataset(se_file_path)

    # Calculate external circulation precipitation
    external_tp = nw_tp + ne_tp + sw_tp + se_tp

    external_tp.to_netcdf(f'/data/external_prev/external_tp/tp_{year}.nc')