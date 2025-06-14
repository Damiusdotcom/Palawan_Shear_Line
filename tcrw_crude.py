import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# Input/output
input_dir = r'Y:\wfs_shared\Personnel Files\DJEV\Palawan_Shear_Line\event'
output_dir = 'tcrw_crude_output'
os.makedirs(output_dir, exist_ok=True)

# Load dataset
ds_tcrw = xr.open_dataset(os.path.join(input_dir, 'tcrw.nc'))
tcrw = ds_tcrw['tcrw']

# Coordinates and times
lon = tcrw.longitude.values
lat = tcrw.latitude.values
times = tcrw.valid_time.values

# Desired extent: [min_lon, max_lon, min_lat, max_lat]
extent = [100, 150, 0, 60]

for i, time in enumerate(times):
    data = tcrw.sel(valid_time=time).values.squeeze()

    fig, ax = plt.subplots(figsize=(10, 6), subplot_kw={'projection': ccrs.PlateCarree()})
    ax.set_extent(extent, crs=ccrs.PlateCarree())

    # Add coastlines and land
    ax.coastlines(resolution='10m')
    ax.add_feature(cfeature.LAND.with_scale('10m'), facecolor='lightgray')
    ax.add_feature(cfeature.BORDERS.with_scale('10m'), linestyle=':')

    # Plot data using pcolormesh for correct geospatial referencing
    # Create 2D lon/lat grids if needed
    lon2d, lat2d = np.meshgrid(lon, lat)
    pcm = ax.pcolormesh(
        lon2d, lat2d, data,
        cmap='Reds',
        vmin=0,
        vmax=2,
        shading='auto',
        transform=ccrs.PlateCarree()
    )

    plt.colorbar(pcm, ax=ax, label='TCRW (kg/m² ≈ mm)', orientation='vertical', fraction=0.046, pad=0.04)
    ax.set_title(f"TCRW - {np.datetime_as_string(time, unit='h')}")

    # Save figure
    filename = f"tcrw_cartopy_{i:03}.png"
    plt.savefig(os.path.join(output_dir, filename), dpi=150, bbox_inches='tight')
    plt.close()
