import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# File paths
mslp_file = r'Y:\wfs_shared\Personnel Files\DJEV\Palawan_Shear_Line\event\mslp.nc'
u10_file = r'Y:\wfs_shared\Personnel Files\DJEV\Palawan_Shear_Line\event\u10.nc'
v10_file = r'Y:\wfs_shared\Personnel Files\DJEV\Palawan_Shear_Line\event\v10.nc'

# Output directory
output_dir = 'mslp_output'
os.makedirs(output_dir, exist_ok=True)

# Load data
mslp_ds = xr.open_dataset(mslp_file)
u10_ds = xr.open_dataset(u10_file)
v10_ds = xr.open_dataset(v10_file)

# Extract variables and convert pressure from Pa to hPa
mslp = mslp_ds['msl'] / 100.0
u10 = u10_ds['u10']
v10 = v10_ds['v10']
times = mslp['valid_time'].values

# Latitude and longitude
lat = mslp['latitude']
lon = mslp['longitude']
lon2d, lat2d = np.meshgrid(lon, lat)

extent = [100, 150, 0, 60]

# Plotting loop
for i, t in enumerate(times):
    mslp_map = mslp.sel(valid_time=t).squeeze()
    u_map = u10.sel(valid_time=t).squeeze()
    v_map = v10.sel(valid_time=t).squeeze()

    fig, ax = plt.subplots(figsize=(10, 8), subplot_kw={'projection': ccrs.PlateCarree()})
    
    # Maintain domain setting
    ax.set_extent(extent, crs=ccrs.PlateCarree())
    
    ax.coastlines(resolution='10m')
    ax.add_feature(cfeature.BORDERS, linewidth=0.5)
    ax.add_feature(cfeature.LAND, facecolor='lightgray')
    ax.add_feature(cfeature.OCEAN, facecolor='lightblue')

    # Plot MSLP contours
    cs = ax.contourf(lon2d, lat2d, mslp_map, levels=np.arange(970, 1040, 1), cmap='coolwarm_r', extend='both')
    cbar = plt.colorbar(cs, ax=ax, orientation='vertical', pad=0.05)
    cbar.set_label('Mean Sea Level Pressure (hPa)')

    # Wind quiver (every nth point)
    skip = 5
    ax.quiver(lon2d[::skip, ::skip], lat2d[::skip, ::skip],
              u_map[::skip, ::skip], v_map[::skip, ::skip],
              scale=700, width=0.0025, color='k')

    # Title
    timestamp = np.datetime_as_string(t, unit='h')
    ax.set_title(f'MSLP and 10 m Wind - {timestamp} UTC')

    # Save figure
    plt.savefig(os.path.join(output_dir, f'mslp_wind_{i:03d}.png'), dpi=150, bbox_inches='tight')
    plt.close()
