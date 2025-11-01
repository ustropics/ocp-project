#!/usr/bin/env python3
# --------------------------------------------------------------
#  Plot & SAVE SSH map with timestamp
#  Robust _FillValue handling + auto-create output dir
# --------------------------------------------------------------

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from datetime import datetime
import os
from netCDF4 import Dataset  # <-- Critical for robust _FillValue

# ------------------------------------------------------------------
# 1. Variables and file paths
# ------------------------------------------------------------------
day = '001'
cmap = 'Spectral_r'

file_path = f'output/data/25deg/2d/gomb4_daily_2024_{day}_2d.nc'
save_dir = 'output/images'
os.makedirs(save_dir, exist_ok=True)
save_path = os.path.join(save_dir, f'ssh_2024_{day}.png')

cmap = 'Spectral_r'

# ------------------------------------------------------------------
# 2. Step 1: Open with netCDF4 to get _FillValue reliably
# ------------------------------------------------------------------
nc = Dataset(file_path, 'r')
fill_value = nc.variables['ssh']._FillValue  # This works even if xarray fails
nc.close()
print(f"_FillValue from netCDF4: {fill_value}")

# ------------------------------------------------------------------
# 3. Step 2: Open with xarray for easy access
# ------------------------------------------------------------------
ds = xr.open_dataset(file_path)
ssh = ds['ssh'].squeeze()
lat = ds['Latitude'].values
lon = ds['Longitude'].values

# ------------------------------------------------------------------
# 4. Human-readable timestamp from 'Date'
# ------------------------------------------------------------------
date_val = ds['Date'].values.item()
date_str = f"{int(date_val):08d}"
date_obj = datetime.strptime(date_str, "%Y%m%d")

month_name = date_obj.strftime("%B")
day = date_obj.day
ordinal = lambda n: "%d%s" % (n, "tsnrhtdd"[(n//10%10!=1)*(n%10<4)*n%10::4])
day_str = ordinal(day)
timestamp = f"{month_name} {day_str}, {date_obj.year}"
print(f"Data date: {timestamp}")

# ------------------------------------------------------------------
# 5. Mask fill value and compute stats
# ------------------------------------------------------------------
ssh_valid = ssh.where(ssh != fill_value)

ssh_min  = float(ssh_valid.min())
ssh_max  = float(ssh_valid.max())
ssh_mean = float(ssh_valid.mean())

print(f"SSH  min : {ssh_min: .4f} m")
print(f"SSH  max : {ssh_max: .4f} m")
print(f"SSH mean : {ssh_mean: .4f} m")

# ------------------------------------------------------------------
# 6. Plot (only the part that changed)
# ------------------------------------------------------------------
fig = plt.figure(figsize=(11, 8))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([lon.min(), lon.max(), lat.min(), lat.max()], crs=ccrs.PlateCarree())

levels = np.linspace(-1.0, 0.5, 41)

# Filled contours
cf = ax.contourf(lon, lat, ssh, levels=levels,
                 cmap=cmap, extend='both',
                 transform=ccrs.PlateCarree())

# *** CONTOUR LINES + LABELS ***
cs = ax.contour(lon, lat, ssh, levels=levels[::5],
                colors='k', linewidths=0.5,
                transform=ccrs.PlateCarree())

# <-- NEW LINE: add the labels
ax.clabel(cs, inline=True, fontsize=9, fmt='%.2f',
          colors='k', inline_spacing=4)

# colour-bar, coastlines, etc. (unchanged)
cbar = plt.colorbar(cf, ax=ax, shrink=0.7, pad=0.06)
cbar.set_label('Sea-Surface Height (m)', fontsize=12)

ax.coastlines(resolution='10m')
ax.add_feature(cfeature.LAND, facecolor='#dddddd')

gl = ax.gridlines(draw_labels=True, linestyle='--', color='gray', alpha=0.5)
gl.top_labels = gl.right_labels = False

ax.set_title(f'Sea-Surface Height ({timestamp})', fontsize=14, pad=15)
ax.text(0.01, 0.01, timestamp, transform=ax.transAxes,
        fontsize=11, color='black',
        bbox=dict(facecolor='white', alpha=0.8, edgecolor='none', pad=3))

# ------------------------------------------------------------------
# 7. Save (no show)
# ------------------------------------------------------------------
plt.tight_layout()
plt.savefig(save_path, dpi=300, bbox_inches='tight')
plt.close(fig)

print(f"Plot saved: {save_path}")