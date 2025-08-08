import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.colors as mcolors
import pandas as pd  # for datetime formatting

# === CONFIGURATION ===
event_tcwv_file = r'event\tcwv.nc'
event_u10_file = r'event\u10.nc'
event_v10_file = r'event\v10.nc'

clim_tcwv_file = r'mean_tcwv_output\mean_tcwv_mean.nc'
clim_wind_file = r'mean_wind_30years_output\30yr_NDJFM_mean_wind.nc'

output_dir = 'tcwv_anomaly_output'
os.makedirs(output_dir, exist_ok=True)

extent = [100, 150, 0, 30]

# Load datasets
ds_tcwv = xr.open_dataset(event_tcwv_file)
ds_u10 = xr.open_dataset(event_u10_file)
ds_v10 = xr.open_dataset(event_v10_file)
clim_tcwv = xr.open_dataset(clim_tcwv_file)['tcwv']
ds_clim_wind = xr.open_dataset(clim_wind_file)
clim_u10 = ds_clim_wind['u10_mean']
clim_v10 = ds_clim_wind['v10_mean']

tcwv = ds_tcwv['tcwv']
u10 = ds_u10['u10']
v10 = ds_v10['v10']

lon = tcwv.longitude.values
lat = tcwv.latitude.values
times = tcwv.valid_time.values

def plot_tcwv_wind_anomaly(tcwv_anom, u10_anom, v10_anom, lon, lat, time, day_index, frame_index, extent, output_dir):
    """Plot TCWV anomaly with wind anomaly quivers and formatted title."""
    lon2d, lat2d = np.meshgrid(lon, lat)

    fig, ax = plt.subplots(figsize=(10, 6), subplot_kw={'projection': ccrs.PlateCarree()})
    ax.set_extent(extent, crs=ccrs.PlateCarree())
    ax.coastlines(resolution='10m')
    ax.add_feature(cfeature.LAND.with_scale('10m'), facecolor='lightgray')
    ax.add_feature(cfeature.BORDERS.with_scale('10m'), linestyle=':')

    lon_min, lon_max, lat_min, lat_max = extent

    # Compute ticks
    lon_ticks = np.arange(np.ceil(lon_min / 10) * 10, lon_max + 1, 10)
    lat_ticks = np.arange(np.ceil(lat_min / 10) * 10, lat_max + 1, 10)

    ax.set_xticks(lon_ticks, crs=ccrs.PlateCarree())
    ax.set_yticks(lat_ticks, crs=ccrs.PlateCarree())

    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, pos: f"{int(v)}°E"))
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, pos: f"{int(v)}°N"))
    ax.tick_params(axis='both', labelsize=10)

    # TCWV anomaly shading
    colors = [
        "#5D3A00", "#B8860B", "#F5DEB3", "#E6E6FA", "#9370DB", "#4B0082"
    ]
    cmap = mcolors.LinearSegmentedColormap.from_list("brown_violet", colors)

    pcm = ax.pcolormesh(
        lon2d, lat2d, tcwv_anom,
        cmap=cmap,
        vmin=-20, vmax=20,
        shading='auto',
        transform=ccrs.PlateCarree()
    )
    plt.colorbar(pcm, ax=ax, label='TCWV Anomaly (mm)', orientation='vertical', fraction=0.046, pad=0.04)

    # Wind anomaly
    skip = 5
    ax.quiver(
        lon[::skip], lat[::skip],
        u10_anom[::skip, ::skip], v10_anom[::skip, ::skip],
        scale=200, width=0.005, color='black',
        transform=ccrs.PlateCarree()
    )

    # Format date and title
    date_str = pd.to_datetime(str(time)).strftime("%d %B %Y")
    day_str = f"+{day_index}" if day_index > 0 else f"{day_index}"
    ax.set_title(f"Day {day_str} ({date_str})", fontsize=14)

    # Save figure (filenames start from 000)
    filename = f"tcwv_wind_anomaly_{frame_index:03}.png"
    plt.savefig(os.path.join(output_dir, filename), dpi=150, bbox_inches='tight')
    plt.close()

# === Loop through time steps ===
for i, time in enumerate(times):
    tcwv_evt = tcwv.sel(valid_time=time).values.squeeze()
    u10_evt = u10.sel(valid_time=time).values
    v10_evt = v10.sel(valid_time=time).values

    tcwv_anom = tcwv_evt - clim_tcwv.values
    u10_anom = u10_evt - clim_u10.values
    v10_anom = v10_evt - clim_v10.values

    day_index = i - 5
    plot_tcwv_wind_anomaly(tcwv_anom, u10_anom, v10_anom, lon, lat, time, day_index, i, extent, output_dir)
