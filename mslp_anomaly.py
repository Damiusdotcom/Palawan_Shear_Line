import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

def plot_mslp_wind_anomaly(mslp_anom, u_anom, v_anom, lon, lat, time, index, extent, output_dir):
    """Plot MSLP anomaly with wind anomaly quivers."""
    lon2d, lat2d = np.meshgrid(lon, lat)

    fig, ax = plt.subplots(figsize=(10, 6), subplot_kw={'projection': ccrs.PlateCarree()})
    ax.set_extent(extent, crs=ccrs.PlateCarree())
    ax.coastlines(resolution='10m')
    ax.add_feature(cfeature.LAND.with_scale('10m'), facecolor='lightgray')
    ax.add_feature(cfeature.BORDERS.with_scale('10m'), linestyle=':')

    # MSLP anomaly shading
    pcm = ax.pcolormesh(
        lon2d, lat2d, mslp_anom,
        cmap='coolwarm_r',
        vmin=-30, vmax=30,
        shading='auto',
        transform=ccrs.PlateCarree()
    )
    plt.colorbar(pcm, ax=ax, label='MSLP Anomaly (hPa)', orientation='vertical', fraction=0.046, pad=0.04)

    # Wind anomaly quiver
    skip = 5
    ax.quiver(
        lon[::skip], lat[::skip],
        u_anom[::skip, ::skip], v_anom[::skip, ::skip],
        scale=200, width=0.005, color='black',
        transform=ccrs.PlateCarree()
    )

    timestamp = np.datetime_as_string(time, unit='h')
    ax.set_title(f"MSLP Anomaly with Wind Anomaly - {timestamp} UTC")

    plt.savefig(os.path.join(output_dir, f"mslp_wind_anomaly_{index:03}.png"), dpi=150, bbox_inches='tight')
    plt.close()

# === Main script ===

# Input/output
input_dir = r'Y:\wfs_shared\Personnel Files\DJEV\Palawan_Shear_Line\event'
clim_mslp_file = r'mean_mslp_output\ndjfm_mslp_mean.nc'
clim_wind_file = r'mean_wind_30years_output\30yr_NDJFM_mean_wind.nc'
output_dir = 'mslp_anomaly_output'
os.makedirs(output_dir, exist_ok=True)

# Load event data
ds_mslp = xr.open_dataset(os.path.join(input_dir, 'mslp.nc'))
ds_u10 = xr.open_dataset(os.path.join(input_dir, 'u10.nc'))
ds_v10 = xr.open_dataset(os.path.join(input_dir, 'v10.nc'))

# Load climatology
clim_mslp = xr.open_dataset(clim_mslp_file)['msl']
clim_u10 = xr.open_dataset(clim_wind_file)['u10_mean']
clim_v10 = xr.open_dataset(clim_wind_file)['v10_mean']

# Variables
mslp = ds_mslp['msl']
u10 = ds_u10['u10']
v10 = ds_v10['v10']

lon = mslp.longitude.values
lat = mslp.latitude.values
times = mslp.valid_time.values

extent = [100, 150, 0, 60]

# Loop through time steps and call the function
for i, time in enumerate(times):
    mslp_evt = mslp.sel(valid_time=time).values.squeeze() / 100.0  # Pa to hPa
    u10_evt = u10.sel(valid_time=time).values
    v10_evt = v10.sel(valid_time=time).values

    mslp_anom = mslp_evt - clim_mslp.values.squeeze() / 100.0
    u_anom = u10_evt - clim_u10.values.squeeze()
    v_anom = v10_evt - clim_v10.values.squeeze()

    plot_mslp_wind_anomaly(mslp_anom, u_anom, v_anom, lon, lat, time, i, extent, output_dir)
