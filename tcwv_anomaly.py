import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.colors as mcolors

# Input/output file paths
event_tcwv_file = r'Y:\wfs_shared\Personnel Files\DJEV\Palawan_Shear_Line\event\tcwv.nc'
event_u10_file = r'Y:\wfs_shared\Personnel Files\DJEV\Palawan_Shear_Line\event\u10.nc'
event_v10_file = r'Y:\wfs_shared\Personnel Files\DJEV\Palawan_Shear_Line\event\v10.nc'

clim_tcwv_file = r'mean_tcwv_output\mean_tcwv_mean.nc'
clim_wind_file = r'mean_wind_30years_output\30yr_NDJFM_mean_wind.nc'

output_dir = 'tcwv_anomaly_output'
os.makedirs(output_dir, exist_ok=True)

# Load event datasets
ds_tcwv = xr.open_dataset(event_tcwv_file)
ds_u10 = xr.open_dataset(event_u10_file)
ds_v10 = xr.open_dataset(event_v10_file)

# Load climatology datasets
clim_tcwv = xr.open_dataset(clim_tcwv_file)['tcwv']
ds_clim_wind = xr.open_dataset(clim_wind_file)
clim_u10 = ds_clim_wind['u10_mean']
clim_v10 = ds_clim_wind['v10_mean']

# Extract event data variables
tcwv = ds_tcwv['tcwv']
u10 = ds_u10['u10']
v10 = ds_v10['v10']

lon = tcwv.longitude.values
lat = tcwv.latitude.values
times = tcwv.valid_time.values

extent = [100, 150, 0, 60]

for i, time in enumerate(times):
    # Event data at this time
    tcwv_evt = tcwv.sel(valid_time=time).values.squeeze()
    u10_evt = u10.sel(valid_time=time).values
    v10_evt = v10.sel(valid_time=time).values

    # Compute anomalies
    tcwv_anom = tcwv_evt - clim_tcwv.values
    u10_anom = u10_evt - clim_u10.values
    v10_anom = v10_evt - clim_v10.values

    # Set up figure and map
    fig, ax = plt.subplots(figsize=(10, 6), subplot_kw={'projection': ccrs.PlateCarree()})
    ax.set_extent(extent, crs=ccrs.PlateCarree())
    ax.coastlines(resolution='10m')
    ax.add_feature(cfeature.LAND.with_scale('10m'), facecolor='lightgray')
    ax.add_feature(cfeature.BORDERS.with_scale('10m'), linestyle=':')

    lon_min, lon_max, lat_min, lat_max = extent

    # Compute ticks strictly within extent, multiples of 5
    lon_ticks = np.arange(np.ceil(lon_min / 10) * 10, lon_max + 1, 10)
    lat_ticks = np.arange(np.ceil(lat_min / 10) * 10, lat_max + 1, 10)

    ax.set_xticks(lon_ticks, crs=ccrs.PlateCarree())
    ax.set_yticks(lat_ticks, crs=ccrs.PlateCarree())

    ax.xaxis.set_major_formatter(plt.FuncFormatter(
        lambda v, pos: f"{int(v)}°E"
    ))
    ax.yaxis.set_major_formatter(plt.FuncFormatter(
        lambda v, pos: f"{int(v)}°N"
    ))

    ax.tick_params(axis='both', labelsize=10)

    # Plot TCWV anomaly
    lon2d, lat2d = np.meshgrid(lon, lat)
    # Define brown to violet colors
    colors = [
        "#5D3A00",  # dark brown (negative extreme)
        "#B8860B",  # lighter brown
        "#F5DEB3",  # pale brown (near zero)
        "#E6E6FA",  # pale violet (near zero)
        "#9370DB",  # medium violet
        "#4B0082"   # dark violet (positive extreme)
    ]

    # Create a ListedColormap and a BoundaryNorm to center zero
    cmap = mcolors.LinearSegmentedColormap.from_list("brown_violet", colors)

    pcm = ax.pcolormesh(
        lon2d, lat2d, tcwv_anom,
        cmap=cmap,
        vmin=-20, vmax=20,
        shading='auto',
        transform=ccrs.PlateCarree()
    )
    plt.colorbar(pcm, ax=ax, label='TCWV Anomaly (mm)', orientation='vertical', fraction=0.046, pad=0.04)

    # Wind anomaly quivers
    skip = 5
    ax.quiver(
        lon[::skip], lat[::skip],
        u10_anom[::skip, ::skip], v10_anom[::skip, ::skip],
        scale=200, width=0.005, color='black',
        transform=ccrs.PlateCarree()
    )

    # Title and save
    timestamp = np.datetime_as_string(time, unit='h')
    ax.set_title(f"TCWV Anomaly with Wind Anomaly - {timestamp} UTC")

    filename = f"tcwv_wind_anomaly_{i:03}.png"
    plt.savefig(os.path.join(output_dir, filename), dpi=150, bbox_inches='tight')
    plt.close()
