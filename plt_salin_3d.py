#!/usr/bin/env python3
"""
plt_salin_3d_stacked_final.py
3D stacked |S × V| at 0, 400, 800 m
SURFACE AT TOP — VERTICAL SPACING WORKS
All layers clearly visible
"""

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from datetime import datetime
import os
from netCDF4 import Dataset

# ------------------------------------------------------------------
# 1. Configuration
# ------------------------------------------------------------------
day = '001'
file_path = f'output/data/25deg/3d/gomb4_daily_2024_{day}_3z.nc'
save_dir = 'output/images'
os.makedirs(save_dir, exist_ok=True)
save_path = os.path.join(
    save_dir, f'salinity_transport_3d_final_2024_{day}.png')

# === 3 DEPTH LEVELS ONLY ===
depth_levels = [0, 50, 400, 800]          # meters
vertical_spacing = 500                    # ← LARGE GAP (controls stretch!)
cmap = 'viridis'

# Domain
lat_min, lat_max = 22.0, 28.0
lon_west, lon_east = -96.0, -84.0

# ------------------------------------------------------------------
# 2. Load valid_range
# ------------------------------------------------------------------
nc = Dataset(file_path, 'r')
valid_s = nc.variables['salinity'].valid_range
valid_u = nc.variables['u'].valid_range
valid_v = nc.variables['v'].valid_range
nc.close()

# ------------------------------------------------------------------
# 3. Load data
# ------------------------------------------------------------------
ds = xr.open_dataset(file_path)
ds = ds.sel(
    Latitude=slice(lat_min, lat_max),
    Longitude=slice(lon_west, lon_east),
    Depth=slice(0, 1000)  # safe upper bound
)

S = ds['salinity'].squeeze()
u = ds['u'].squeeze()
v = ds['v'].squeeze()
lat = ds['Latitude'].values
lon = ds['Longitude'].values
depth_full = ds['Depth'].values

# ------------------------------------------------------------------
# 4. Timestamp
# ------------------------------------------------------------------
date_val = ds['Date'].values.item()
date_str = f"{int(date_val):08d}"
date_obj = datetime.strptime(date_str, "%Y%m%d")


def ordinal(n):
    return "%d%s" % (n, "tsnrhtdd"[(n//10 % 10 != 1)*(n % 10 < 4)*n % 10::4])


timestamp = f"{date_obj.strftime('%B')} {ordinal(date_obj.day)}, {date_obj.year}"

# ------------------------------------------------------------------
# 5. Mask & compute |S × V|
# ------------------------------------------------------------------
S = S.where((S >= valid_s[0]) & (S <= valid_s[1]), np.nan)
u = u.where((u >= valid_u[0]) & (u <= valid_u[1]), np.nan)
v = v.where((v >= valid_v[0]) & (v <= valid_v[1]), np.nan)

transport_mag = np.sqrt((S * u)**2 + (S * v)**2)  # (Depth, Lat, Lon)

# ------------------------------------------------------------------
# 6. Prepare 3D plot
# ------------------------------------------------------------------
fig = plt.figure(figsize=(25, 25))  # tall figure
ax = fig.add_subplot(111, projection='3d')

Lon2D, Lat2D = np.meshgrid(lon, lat)

# Color scaling
all_vals = []
for d in depth_levels:
    if d in depth_full:
        idx = np.abs(depth_full - d).argmin()
        slice_data = transport_mag.isel(Depth=idx).values
        all_vals.extend(slice_data[~np.isnan(slice_data)])
vmin = 0
vmax = np.percentile(all_vals, 98) if all_vals else 40

# ------------------------------------------------------------------
# 7. PLOT: 3 LEVELS — VERTICAL SPACING WORKS
# ------------------------------------------------------------------
base_z = len(depth_levels) * vertical_spacing  # e.g. 2400

for i, target_depth in enumerate(depth_levels):
    if target_depth not in depth_full:
        print(f"Warning: {target_depth} m not in data. Skipping.")
        continue

    idx = np.abs(depth_full - target_depth).argmin()
    Z_slice = transport_mag.isel(Depth=idx).values  # (Lat, Lon)

    z_level = base_z - i * vertical_spacing  # 2400, 1600, 800

    norm_data = np.clip((Z_slice - vmin) / (vmax - vmin), 0, 1)
    alpha = 0.95 if i == 0 else 0.85  # surface = solid

    ax.plot_surface(
        Lon2D, Lat2D, np.full_like(Z_slice, z_level),
        facecolors=plt.cm.viridis(norm_data),
        rstride=1, cstride=1,
        linewidth=0.4, edgecolor='k', alpha=alpha, shade=False
    )

# ------------------------------------------------------------------
# 8. Z-AXIS: Show actual depth (0 m at top)
# ------------------------------------------------------------------
z_ticks = [base_z - i * vertical_spacing for i in range(len(depth_levels))]
ax.set_zticks(z_ticks)
ax.set_zticklabels([f"{d} m" for d in depth_levels])
ax.set_zlabel('Depth (m)', fontsize=12, labelpad=15)

# CRITICAL: Z-LIMIT MUST MATCH SPACING!
ax.set_zlim(0, base_z + 200)

# ------------------------------------------------------------------
# 9. Axes labels & view
# ------------------------------------------------------------------
ax.set_xlabel('Longitude (°W)', fontsize=12, labelpad=15)
ax.set_ylabel('Latitude (°N)', fontsize=12, labelpad=15)
ax.set_xlim(lon.min(), lon.max())
ax.set_ylim(lat.min(), lat.max())

# BEST VIEW: See all layers
ax.view_init(elev=15, azim=-60)

# Title
ax.set_title(
    f'3D Stacked Salinity Transport |S × V|\n'
    f'0 m (top) → 800 m (bottom) | {timestamp}',
    fontsize=14, pad=30
)

# ------------------------------------------------------------------
# 10. COLORBAR — TIGHTER TO FIGURE
# ------------------------------------------------------------------
m = plt.cm.ScalarMappable(cmap=cmap)
m.set_array([vmin, vmax])
m.set_clim(vmin, vmax)

# TIGHTER COLORBAR: shrink + less pad
cbar = fig.colorbar(m, ax=ax, shrink=0.5, pad=0.03, aspect=20)
cbar.set_label('Salinity Transport |S × V| (psu·m/s)', fontsize=11)

# ------------------------------------------------------------------
# TITLE — LESS WHITE SPACE
# ------------------------------------------------------------------
ax.set_title(
    f'3D Stacked Salinity Transport |S × V|\n'
    f'0 m (top) → 800 m (bottom) | {timestamp}',
    fontsize=16, pad=10  # ← reduced from 30
)

# ------------------------------------------------------------------
# 11. Save
# ------------------------------------------------------------------
plt.tight_layout()
plt.savefig(save_path, dpi=300, bbox_inches='tight')
plt.close(fig)
print(f"SUCCESS: 3D plot saved with full vertical stretch: {save_path}")
