#!/usr/bin/env python3
"""
plt_salin_eddy_combined_no_flow.py
COMBINES 2D + 3D datasets:
  - 3D: salinity + velocity → |S × V|
  - 2D: ssh → eddy detection
NO BAROTROPIC FLOW (quiver removed)
"""

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from datetime import datetime
import os
from netCDF4 import Dataset

# ------------------------------------------------------------------
# 1. Configuration
# ------------------------------------------------------------------
day = '001'
cmap = 'cubehelix'

# FILE PATHS
file_3d = f'output/data/25deg/3d/gomb4_daily_2024_{day}_3z.nc'
# ← adjust if needed
file_2d = f'output/data/25deg/2d/gomb4_daily_2024_{day}_2d.nc'

save_dir = 'output/images'
os.makedirs(save_dir, exist_ok=True)
save_path = os.path.join(
    save_dir, f'salinity_transport_eddy_combined_2024_{day}.png')

# ------------------------------------------------------------------
# 2. Load 3D data (salinity + velocity)
# ------------------------------------------------------------------
nc3 = Dataset(file_3d, 'r')
valid_s = nc3.variables['salinity'].valid_range
valid_u = nc3.variables['u'].valid_range
valid_v = nc3.variables['v'].valid_range
nc3.close()

ds3 = xr.open_dataset(file_3d)
S = ds3['salinity'].squeeze()
u = ds3['u'].squeeze()
v = ds3['v'].squeeze()
lat = ds3['Latitude'].values
lon = ds3['Longitude'].values

# Mask
S = S.where((S >= valid_s[0]) & (S <= valid_s[1]), np.nan)
u = u.where((u >= valid_u[0]) & (u <= valid_u[1]), np.nan)
v = v.where((v >= valid_v[0]) & (v <= valid_v[1]), np.nan)

# Depth-averaged |S × V|
Su = S * u
Sv = S * v
Su_mean = Su.mean(dim='Depth', skipna=True)
Sv_mean = Sv.mean(dim='Depth', skipna=True)
transport_mag = np.sqrt(Su_mean**2 + Sv_mean**2)
transport_mag = transport_mag.where(np.isfinite(transport_mag), np.nan)

mag_max = float(transport_mag.max()) if transport_mag.notnull().any() else 40

# ------------------------------------------------------------------
# 3. Load 2D data (SSH only)
# ------------------------------------------------------------------
ds2 = xr.open_dataset(file_2d)
ssh = ds2['ssh'].squeeze()

# Ensure same grid
assert np.allclose(ds2['Latitude'].values, lat), "Latitude mismatch!"
assert np.allclose(ds2['Longitude'].values, lon), "Longitude mismatch!"

# ------------------------------------------------------------------
# 4. Timestamp
# ------------------------------------------------------------------
date_val = ds3['Date'].values.item()
date_str = f"{int(date_val):08d}"
date_obj = datetime.strptime(date_str, "%Y%m%d")


def ordinal(n):
    return "%d%s" % (n, "tsnrhtdd"[(n//10 % 10 != 1)*(n % 10 < 4)*n % 10::4])


timestamp = f"{date_obj.strftime('%B')} {ordinal(date_obj.day)}, {date_obj.year}"

# ------------------------------------------------------------------
# 5. Plot
# ------------------------------------------------------------------
fig = plt.figure(figsize=(13, 10))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([lon.min(), lon.max(), lat.min(),
              lat.max()], crs=ccrs.PlateCarree())

# 1. Color: |S × V| from 3D
levels = np.linspace(0, 50, 51)
cf = ax.contourf(lon, lat, transport_mag, levels=levels,
                 cmap=cmap, extend='max', transform=ccrs.PlateCarree())

# 2. Contours: SSH from 2D (eddy detection)
ssh_levels = np.linspace(-0.6, 0.4, 13)
ssh_levels2 = np.linspace(.4, 0.6, 7)

cs = ax.contour(lon, lat, ssh, levels=ssh_levels,
                colors='white', linewidths=0.9, alpha=0.9,
                transform=ccrs.PlateCarree())
ax.clabel(cs, inline=True, fontsize=8, fmt='%.2f', colors='white')

# NO QUIVER — BAROTROPIC FLOW REMOVED

# Colorbar
cbar = plt.colorbar(cf, ax=ax, shrink=0.7, pad=0.06)
cbar.set_label('Depth-Averaged |S × V| (psu·m/s)', fontsize=11)

# Map features
ax.coastlines(resolution='10m', linewidth=1, color='black')
ax.add_feature(cfeature.LAND, facecolor="#d4a574")
ax.add_feature(cfeature.BORDERS, linewidth=0.5, color='black')
ax.add_feature(cfeature.LAKES, facecolor='#a0c8f0',
               edgecolor='black', linewidth=0.5)
ax.add_feature(cfeature.STATES, linewidth=0.5, edgecolor='gray')

# Grid
gl = ax.gridlines(draw_labels=True, linestyle='--', color='gray', alpha=0.5)
gl.top_labels = gl.right_labels = False

# Title
ax.set_title(
    f'Eddy-Driven Salinity Transport ({timestamp})\n'
    f'Depth-Averaged Salinity Transport (shaded, |S × V|) | Sea Surface Heights (contoured, meters)',
    fontsize=13, pad=15
)

# Stats
ax.text(0.01, 0.01,
        f"Max |S×V|: {mag_max:.3f} psu·m/s",
        transform=ax.transAxes, fontsize=10, color='white',
        bbox=dict(facecolor='black', alpha=0.8, edgecolor='none', pad=3))

# ------------------------------------------------------------------
# 6. Save
# ------------------------------------------------------------------
plt.tight_layout()
plt.savefig(save_path, dpi=300, bbox_inches='tight')
plt.close(fig)
print(f"Plot saved (no barotropic flow): {save_path}")
