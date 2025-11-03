#!/usr/bin/env python3
"""
plt_salin_cross_section.py
Standalone script: Latitudinal cross-section of depth-averaged salinity transport
Averaged over longitude band: 90°W → 87°W
Saves: output/images/salinity_transport_cross_2024_XXX.png
"""

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import os
from netCDF4 import Dataset

# ------------------------------------------------------------------
# 1. Configuration
# ------------------------------------------------------------------
cmap = 'cubehelix'

file_path = f'output/data/25deg/3d/gomb4_daily_2024_{day}_3z.nc'
save_dir = 'output/images'
os.makedirs(save_dir, exist_ok=True)
save_path = os.path.join(save_dir, f'salinity_transport_cross_2024_{day}.png')

# Longitude band for cross-section (90°W → 87°W)
lon_west, lon_east = -90.0, -87.0

# ------------------------------------------------------------------
# 2. Load valid_range from NetCDF (robust masking)
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
# 3. Load data with xarray
# ------------------------------------------------------------------
ds = xr.open_dataset(file_path)

# Squeeze time dimension
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


def ordinal(n):
    return "%d%s" % (n, "tsnrhtdd"[(n//10 % 10 != 1)*(n % 10 < 4)*n % 10::4])


timestamp = f"{month_name} {ordinal(date_obj.day)}, {date_obj.year}"
print(f"Data date: {timestamp}")

# ------------------------------------------------------------------
# 5. Apply VALID_RANGE masking
# ------------------------------------------------------------------
S = S.where((S >= valid_s[0]) & (S <= valid_s[1]), np.nan)
u = u.where((u >= valid_u[0]) & (u <= valid_u[1]), np.nan)
v = v.where((v >= valid_v[0]) & (v <= valid_v[1]), np.nan)

# ------------------------------------------------------------------
# 6. Compute salinity transport |S × V|
# ------------------------------------------------------------------
Su = S * u
Sv = S * v

# Depth-average
Su_mean = Su.mean(dim='Depth', skipna=True)
Sv_mean = Sv.mean(dim='Depth', skipna=True)
transport_mag = np.sqrt(Su_mean**2 + Sv_mean**2)

# Ensure finite values
transport_mag = transport_mag.where(np.isfinite(transport_mag), np.nan)

# ------------------------------------------------------------------
# 7. Select longitude band and average zonally
# ------------------------------------------------------------------
lon_mask = (lon >= lon_west) & (lon <= lon_east)
transport_slice = transport_mag.where(lon_mask)  # Mask outside band
transport_zonal = transport_slice.mean(
    dim='Longitude', skipna=True)  # 1D: (lat)

# Stats
valid_vals = transport_zonal.dropna('Latitude')
if len(valid_vals) > 0:
    t_min = float(valid_vals.min())
    t_max = float(valid_vals.max())
    t_mean = float(valid_vals.mean())
else:
    t_min = t_max = t_mean = np.nan

print(f"Cross-section |S×V|  min : {t_min:.4f} psu·m/s")
print(f"Cross-section |S×V|  max : {t_max:.4f} psu·m/s")
print(f"Cross-section |S×V| mean : {t_mean:.4f} psu·m/s")

# ------------------------------------------------------------------
# 8. Plot: Latitudinal Cross-Section
# ------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(10, 6))

# Dynamic contour levels (0 to max, or cap at 40)
levels = np.linspace(0, max(40, np.nanpercentile(valid_vals, 95)), 41)

# Line plot (clear and publication-ready)
ax.plot(lat, transport_zonal, color='darkmagenta',
        linewidth=2.8, label='Depth-averaged |S × V|')

# Optional: filled area under curve
ax.fill_between(lat, transport_zonal, 0, color='mediumpurple', alpha=0.3)

# Styling
ax.set_xlabel('Latitude (°N)', fontsize=13)
ax.set_ylabel('Salinity Transport Magnitude |S × V| (psu·m/s)', fontsize=13)
ax.set_title(
    f'Latitudinal Cross-Section of Depth-Averaged Salinity Transport\n'
    f'Zonally Averaged {lon_west:.0f}°W → {abs(lon_east):.0f}°W  |  {timestamp}',
    fontsize=14, pad=15
)

ax.grid(True, linestyle='--', alpha=0.6)
ax.set_xlim(lat.min(), lat.max())
ax.legend(fontsize=11)

# Add text box with stats
stats_text = f"Min: {t_min:.3f}\nMean: {t_mean:.3f}\nMax: {t_max:.3f}"
ax.text(0.02, 0.98, stats_text, transform=ax.transAxes,
        fontsize=10, verticalalignment='top',
        bbox=dict(boxstyle="round,pad=0.4", facecolor="white", alpha=0.8))

# ------------------------------------------------------------------
# 9. Save
# ------------------------------------------------------------------
plt.tight_layout()
plt.savefig(save_path, dpi=300, bbox_inches='tight')
plt.close(fig)
print(f"Cross-section plot saved: {save_path}")
