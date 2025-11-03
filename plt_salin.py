#!/usr/bin/env python3
# --------------------------------------------------------------
#  FIXED: Depth-averaged Salinity Transport |S × V|
#  Uses VALID_RANGE masking (robust to ncra averaging)
#  No more inf! Saves: output/images/salinity_transport_2024_001.png
# --------------------------------------------------------------

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from datetime import datetime
import os
from netCDF4 import Dataset

# ------------------------------------------------------------------
# 1. File paths
# ------------------------------------------------------------------
cmap = 'cubehelix'

file_path = f'output/data/25deg/3d/gomb4_daily_2024_{day}_3z.nc'
save_dir = 'output/images'
os.makedirs(save_dir, exist_ok=True)
save_path = os.path.join(save_dir, f'salinity_transport_2024_{day}.png')

# ------------------------------------------------------------------
# 2. Get VALID_RANGE (robust masking!)
# ------------------------------------------------------------------
nc = Dataset(file_path, 'r')
valid_s = nc.variables['salinity'].valid_range
valid_u = nc.variables['u'].valid_range
valid_v = nc.variables['v'].valid_range
nc.close()
print(f"Valid ranges - S: [{valid_s[0]:.4f}, {valid_s[1]:.4f}] psu")
print(f"Valid ranges - u: [{valid_u[0]:.4f}, {valid_u[1]:.4f}] m/s")
print(f"Valid ranges - v: [{valid_v[0]:.4f}, {valid_v[1]:.4f}] m/s")

# ------------------------------------------------------------------
# 3. Load with xarray
# ------------------------------------------------------------------
ds = xr.open_dataset(file_path)

# Squeeze time
S = ds['salinity'].squeeze()   # (Depth, Lat, Lon)
u = ds['u'].squeeze()
v = ds['v'].squeeze()
lat = ds['Latitude'].values
lon = ds['Longitude'].values

# ------------------------------------------------------------------
# 4. Timestamp
# ------------------------------------------------------------------
date_val = ds['Date'].values.item()
date_str = f"{int(date_val):08d}"
date_obj = datetime.strptime(date_str, "%Y%m%d")
month_name = date_obj.strftime("%B")
day = date_obj.day


def ordinal(n): return "%d%s" % (
    n, "tsnrhtdd"[(n//10 % 10 != 1)*(n % 10 < 4)*n % 10::4])


timestamp = f"{month_name} {ordinal(day)}, {date_obj.year}"
print(f"Data date: {timestamp}")

# ------------------------------------------------------------------
# 5. MASK using VALID_RANGE (fixes inf!)
# ------------------------------------------------------------------
S = S.where((S >= valid_s[0]) & (S <= valid_s[1]), np.nan)
u = u.where((u >= valid_u[0]) & (u <= valid_u[1]), np.nan)
v = v.where((v >= valid_v[0]) & (v <= valid_v[1]), np.nan)

# ------------------------------------------------------------------
# 6. Salinity transport
# ------------------------------------------------------------------
Su = S * u  # psu * m/s
Sv = S * v

# Depth-average
Su_mean = Su.mean(dim='Depth', skipna=True)
Sv_mean = Sv.mean(dim='Depth', skipna=True)

# Magnitude |S × V|
transport_mag = np.sqrt(Su_mean**2 + Sv_mean**2)

# Force finite (safety net)
transport_mag = transport_mag.where(np.isfinite(transport_mag), np.nan)

# ------------------------------------------------------------------
# 7. Stats (now finite!)
# ------------------------------------------------------------------
mag_valid = transport_mag.stack(z=('Latitude', 'Longitude')).dropna('z')
if len(mag_valid) > 0:
    mag_min = float(mag_valid.min())
    mag_max = float(mag_valid.max())
    mag_mean = float(mag_valid.mean())
else:
    mag_min = mag_max = mag_mean = np.nan

print(f"Salinity transport |S×V|  min : {mag_min:.4f} psu·m/s")
print(f"Salinity transport |S×V|  max : {mag_max:.4f} psu·m/s")
print(f"Salinity transport |S×V| mean : {mag_mean:.4f} psu·m/s")

# ------------------------------------------------------------------
# 8. Plot (dynamic levels based on max)
# ------------------------------------------------------------------
fig = plt.figure(figsize=(11, 8))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([lon.min(), lon.max(), lat.min(),
              lat.max()], crs=ccrs.PlateCarree())

# Dynamic levels: 0 to 95% of max, 31 levels
levels = np.linspace(0, 40, 41)
levels_white = np.linspace(0, 30, 4)
levels_black = np.linspace(30, 40, 2)

cf = ax.contourf(lon, lat, transport_mag, levels=levels,
                 cmap=cmap, extend='max',
                 transform=ccrs.PlateCarree())

cs = ax.contour(lon, lat, transport_mag, levels_white,
                colors='white', linewidths=0.5,
                transform=ccrs.PlateCarree())

cs2 = ax.contour(lon, lat, transport_mag, levels_black,
                 colors='black', linewidths=0.5,
                 transform=ccrs.PlateCarree())

ax.clabel(cs, inline=True, fontsize=9, fmt='%.2f', colors='white')
ax.clabel(cs2, inline=True, fontsize=9, fmt='%.2f', colors='black')

cbar = plt.colorbar(cf, ax=ax, shrink=0.7, pad=0.06)
cbar.set_label('Salinity Transport Magnitude |S × V| (psu·m/s)', fontsize=12)

ax.coastlines(resolution='10m', linewidth=1, color='black', zorder=3)
ax.add_feature(cfeature.LAND, facecolor="#9c6b00", zorder=2)
ax.add_feature(cfeature.BORDERS, linewidth=0.5, zorder=3, color='black')
ax.add_feature(cfeature.LAKES, facecolor='#a0c8f0',
               edgecolor='black', linewidth=0.5, zorder=3)

# add state borders
ax.add_feature(cfeature.STATES, linewidth=0.5, zorder=3,
               facecolor='none', edgecolor='black')

gl = ax.gridlines(draw_labels=True, linestyle='--', color='gray', alpha=0.5)
gl.top_labels = gl.right_labels = False

ax.set_title(
    f'Depth-Averaged Salinity Transport ({timestamp})', fontsize=14, pad=15)
ax.text(0.01, 0.01, timestamp, transform=ax.transAxes,
        fontsize=11, color='white',
        bbox=dict(facecolor='black', alpha=0.7, edgecolor='none', pad=3))

# ------------------------------------------------------------------
# 9. Save
# ------------------------------------------------------------------
plt.tight_layout()
plt.savefig(save_path, dpi=300, bbox_inches='tight')
plt.close(fig)
print(f"Plot saved: {save_path}")
