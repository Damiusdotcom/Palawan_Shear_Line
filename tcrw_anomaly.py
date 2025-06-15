import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# Input/output file paths
event_tcrw_file = r'Y:\wfs_shared\Personnel Files\DJEV\Palawan_Shear_Line\event\tcrw.nc'
event_u10_file = r'Y:\wfs_shared\Personnel Files\DJEV\Palawan_Shear_Line\event\u10.nc'
event_v10_file = r'Y:\wfs_shared\Personnel Files\DJEV\Palawan_Shear_Line\event\v10.nc'

clim_tcrw_file = r'mean_tcwv_output\mean_tcrw_mean.nc'
clim_wind_file = r'mean_wind_30years_output\30yr_NDJFM_mean_wind.nc'

output_dir = 'tcrw_anomaly_output'
os.makedirs(output_dir, exist_ok=True)

# Load event datasets
ds_tcrw = xr.open_dataset(event_tcrw_file)
ds_u10 = xr.open_dataset(event_u10_file)
ds_v10 = xr.open_dataset(event_v10_file)

# Load climatology datasets
clim_tcrw = xr.open_dataset(clim_tcrw_file)['tcrw']
ds_clim_wind = xr.open_dataset(clim_wind_file)
clim_u10 = ds_clim_wind['u10_mean']
clim_v10 = ds_clim_wind['v10_mean']

# Extract event data variables
tcrw = ds_tcrw['tcrw']
u10 = ds_u10['u10']
v10 = ds_v10['v10']

lon = tcrw.longitude.values
lat = tcrw.latitude.values
times = tcrw.valid_time.values

extent = [100, 150, 0, 60]

for i, time in enumerate(times):
    # Event data at this time
    tcrw_evt = tcrw.sel(valid_time=time).values.squeeze()
    u10_evt = u10.sel(valid_time=time).values
    v10_evt = v10.sel(valid_time=time).values

    # Compute anomalies
    tcrw_anom = tcrw_evt - clim_tcrw.values
    u10_anom = u10_evt - clim_u10.values
    v10_anom = v10_evt - clim_v10.values

    # Set up figure and map
    fig, ax = plt.subplots(figsize=(10, 6), subplot_kw={'projection': ccrs.PlateCarree()})
    ax.set_extent(extent, crs=ccrs.PlateCarree())
    ax.coastlines(resolution='10m')
    ax.add_feature(cfeature.LAND.with_scale('10m'), facecolor='lightgray')
    ax.add_feature(cfeature.BORDERS.with_scale('10m'), linestyle=':')

    # Plot TCRW anomaly
    lon2d, lat2d = np.meshgrid(lon, lat)
    pcm = ax.pcolormesh(
        lon2d, lat2d, tcrw_anom,
        cmap='bwr',
        vmin=-1, vmax=1,
        shading='auto',
        transform=ccrs.PlateCarree()
    )
    plt.colorbar(pcm, ax=ax, label='TCRW Anomaly (kg/m²)', orientation='vertical', fraction=0.046, pad=0.04)

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
    ax.set_title(f"TCRW Anomaly with Wind Anomaly - {timestamp} UTC")

    filename = f"tcrw_wind_anomaly_{i:03}.png"
    plt.savefig(os.path.join(output_dir, filename), dpi=150, bbox_inches='tight')
    plt.close()
