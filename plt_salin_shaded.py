#!/usr/bin/env python3
# --------------------------------------------------------------
#  FIXED: Depth-averaged Salinity Transport |S × V|
#  + 50 % shaded longitude bands
#  + TOP-OF-MAP LABEL BOXES ALIGNED WITH ACTUAL LONGITUDE
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
# 2. Get VALID_RANGE
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
S = ds['salinity'].squeeze()
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


def ordinal(n):
    return "%d%s" % (n, "tsnrhtdd"[(n//10 % 10 != 1)*(n % 10 < 4)*n % 10::4])


timestamp = f"{month_name} {ordinal(date_obj.day)}, {date_obj.year}"

# ------------------------------------------------------------------
# 5. Mask & compute transport
# ------------------------------------------------------------------
S = S.where((S >= valid_s[0]) & (S <= valid_s[1]), np.nan)
u = u.where((u >= valid_u[0]) & (u <= valid_u[1]), np.nan)
v = v.where((v >= valid_v[0]) & (v <= valid_v[1]), np.nan)

Su = S * u
Sv = S * v
transport_mag = np.sqrt(Su**2 + Sv**2).mean(dim='Depth', skipna=True)
transport_mag = transport_mag.where(np.isfinite(transport_mag), np.nan)

# ------------------------------------------------------------------
# 6. Stats
# ------------------------------------------------------------------
mag_valid = transport_mag.stack(z=('Latitude', 'Longitude')).dropna('z')
mag_min = float(mag_valid.min()) if len(mag_valid) > 0 else np.nan
mag_max = float(mag_valid.max()) if len(mag_valid) > 0 else np.nan
mag_mean = float(mag_valid.mean()) if len(mag_valid) > 0 else np.nan

# ------------------------------------------------------------------
# 7. Plot
# ------------------------------------------------------------------
fig = plt.figure(figsize=(11, 8))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([lon.min(), lon.max(), lat.min(),
              lat.max()], crs=ccrs.PlateCarree())

# ------------------------------------------------------------------
# 7.1  SHADE BANDS (50% alpha)
# ------------------------------------------------------------------
bands = [
    {"west": -96.0, "east": -93.0, "label": "96°W–93°W"},
    {"west": -93.0, "east": -90.0, "label": "93°W–90°W"},
    {"west": -90.0, "east": -87.0, "label": "90°W–87°W"},
    {"west": -87.0, "east": -84.0, "label": "87°W–84°W"},
]

band_cmap = plt.cm.get_cmap('tab10', len(bands))

for i, b in enumerate(bands):
    rect = plt.Rectangle(
        (b["west"], lat.min()),
        b["east"] - b["west"],
        lat.max() - lat.min(),
        transform=ccrs.PlateCarree(),
        facecolor=band_cmap(i),
        alpha=0.5,
        zorder=1,
        edgecolor='none'
    )
    ax.add_patch(rect)

# ------------------------------------------------------------------
# 7.2  CONTOUR FILL + CONTOURS
# ------------------------------------------------------------------
levels = np.linspace(0, 40, 41)
cf = ax.contourf(lon, lat, transport_mag, levels=levels,
                 cmap=cmap, extend='max', transform=ccrs.PlateCarree())

cs = ax.contour(lon, lat, transport_mag, levels=np.linspace(0, 30, 4),
                colors='white', linewidths=0.5, transform=ccrs.PlateCarree())
cs2 = ax.contour(lon, lat, transport_mag, levels=np.linspace(30, 40, 2),
                 colors='black', linewidths=0.5, transform=ccrs.PlateCarree())

ax.clabel(cs, inline=True, fontsize=9, fmt='%.2f', colors='white')
ax.clabel(cs2, inline=True, fontsize=9, fmt='%.2f', colors='black')

cbar = plt.colorbar(cf, ax=ax, shrink=0.7, pad=0.06)
cbar.set_label('Salinity Transport Magnitude |S × V| (psu·m/s)', fontsize=12)

# ------------------------------------------------------------------
# 7.3  MAP FEATURES
# ------------------------------------------------------------------
ax.coastlines(resolution='10m', linewidth=1, color='black', zorder=3)
ax.add_feature(cfeature.LAND, facecolor="#9c6b00", zorder=2)
ax.add_feature(cfeature.BORDERS, linewidth=0.5, zorder=3, color='black')
ax.add_feature(cfeature.LAKES, facecolor='#a0c8f0',
               edgecolor='black', linewidth=0.5, zorder=3)
ax.add_feature(cfeature.STATES, linewidth=0.5, zorder=3,
               facecolor='none', edgecolor='black')

gl = ax.gridlines(draw_labels=True, linestyle='--', color='gray', alpha=0.5)
gl.top_labels = gl.right_labels = False

ax.set_title(
    f'Depth-Averaged Salinity Transport ({timestamp})', fontsize=14, pad=50)

ax.text(0.01, 0.01, timestamp, transform=ax.transAxes,
        fontsize=11, color='white',
        bbox=dict(facecolor='black', alpha=0.7, edgecolor='none', pad=3))

# ------------------------------------------------------------------
# 7.4  TOP LABEL BOXES ALIGNED WITH LONGITUDE
# ------------------------------------------------------------------
# Helper: convert (lon, lat) to figure fraction


def lon_to_figx(lon_val):
    x_display = ax.transData.transform((lon_val, lat.min()))[0]
    return fig.transFigure.inverted().transform((x_display, 0))[0]


# Height and vertical position of label boxes
box_height = 0.028
box_y = ax.get_position().y1 + 0.038

for i, b in enumerate(bands):
    x_left = lon_to_figx(b["west"])
    x_right = lon_to_figx(b["east"])
    box_width = x_right - x_left

    label_ax = fig.add_axes(
        [x_left, box_y, box_width, box_height], frameon=False)
    label_ax.set_xlim(0, 1)
    label_ax.set_ylim(0, 1)
    label_ax.axis('off')

    # Colored background
    label_ax.add_patch(
        plt.Rectangle((0.05, 0.2), 0.9, 0.6,
                      facecolor=band_cmap(i),
                      edgecolor='black', linewidth=0.8,
                      transform=label_ax.transAxes)
    )
    # Centered text
    label_ax.text(0.5, 0.5, b["label"],
                  ha='center', va='center',
                  fontsize=9, fontweight='bold',
                  transform=label_ax.transAxes)

# ------------------------------------------------------------------
# 8. Final layout & save
# ------------------------------------------------------------------
plt.subplots_adjust(top=0.88, bottom=0.1, left=0.05, right=0.95)
plt.tight_layout()
plt.savefig(save_path, dpi=300, bbox_inches='tight')
plt.close(fig)
print(f"Plot saved: {save_path}")
