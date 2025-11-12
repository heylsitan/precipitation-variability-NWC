import xarray as xr
import numpy as np
import os
import pandas as pd
from glob import glob

# Date range: All summer days from 1982–2020
date_list = pd.date_range("19820101", "20201231")
summer_dates = [date for date in date_list if date.month in [6, 7, 8]]

# Constants
R = 6.371e6  # Earth radius (meters)
g = 9.8      # Gravitational acceleration (m/s²)

for date in summer_dates:
    # One file per day
    file_pattern = f"/data2/MERRA2/3d/{date.year}/MERRA2_*.inst3_3d_asm_Np.{date.strftime('%Y%m%d')}.nc4"
    file_list = sorted(glob(file_pattern))

    if not file_list:
        print(f"No files found for {date.strftime('%Y-%m-%d')}")
        continue

    # === Read data ===
    ds = xr.open_dataset(file_list[0])

    # Extract variables
    q = ds['QV']   # Specific humidity [time, level, lat, lon]
    u = ds['U']
    v = ds['V']
    lev = ds['lev']
    p = lev * 100  # Pa

    # === Fix level order (ensure from high altitude → low altitude) ===
    if lev.values[0] > lev.values[-1]:  # ✅ Note to add .values
        lev = lev[::-1]
        q = q.sel(lev=lev)
        u = u.sel(lev=lev)
        v = v.sel(lev=lev)
        p = p[::-1]

    # === Calculate dp ===
    dp = np.gradient(p)  # ✅ From top to bottom dp is positive
    dp_da = xr.DataArray(dp, coords={'lev': lev}, dims=['lev'])

    # === Calculate total water vapor flux (vertical integration) ===
    qu_all = q * u
    qv_all = q * v
    qu_int = (qu_all * dp_da).sum(dim='lev') / g
    qv_int = (qv_all * dp_da).sum(dim='lev') / g

    # === Daily time average (optional) ===
    qu = qu_int.mean(dim='time')
    qv = qv_int.mean(dim='time')
    qu.name = 'qu'
    qv.name = 'qv'

    # === Longitude and latitude grid information ===
    lat = ds['lat']
    lon = ds['lon']
    lat_rad = np.deg2rad(lat)
    lon_rad = np.deg2rad(lon)

    dlat = np.gradient(lat_rad)
    dlon = np.gradient(lon_rad)
    dlat_da = xr.DataArray(dlat, coords={'lat': lat}, dims=['lat'])
    dlon_da = xr.DataArray(dlon, coords={'lon': lon}, dims=['lon'])
    cos_lat = xr.DataArray(np.cos(lat_rad), coords={'lat': lat}, dims=['lat'])

    cos_lat_2d, _ = xr.broadcast(cos_lat, lon)
    dlat_2d, _ = xr.broadcast(dlat_da, lon)
    dlon_2d, _ = xr.broadcast(dlon_da, lat)

    dy = R * dlat_2d
    dx = R * cos_lat_2d * dlon_2d

    # === Water vapor flux divergence ===
    dqu_dx = qu.diff('lon') / dx.isel(lon=slice(1, None))
    dqv_dy = qv.diff('lat') / dy.isel(lat=slice(1, None))
    dqu_dx = dqu_dx.pad(lon=(1, 0), constant_values=np.nan)
    dqv_dy = dqv_dy.pad(lat=(1, 0), constant_values=np.nan)
    divergence = dqu_dx + dqv_dy
    divergence.name = 'divQ'

    # === Add time dimension ===
    current_time = pd.Timestamp(date)
    qu = qu.expand_dims(time=[current_time])
    qv = qv.expand_dims(time=[current_time])
    divergence = divergence.expand_dims(time=[current_time])

    # === Merge and save ===
    out_ds = xr.merge([qu, qv, divergence])
    out_path = f"/data1/moisture_flux/moisture_flux_{date.strftime('%Y-%m-%d')}.nc"
    out_ds.to_netcdf(out_path)

    print(f"{date.strftime('%Y-%m-%d')} done.")