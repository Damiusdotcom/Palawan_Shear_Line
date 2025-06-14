import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

input_dir = r'Y:\wfs_shared\Personnel Files\DJEV\Palawan_Shear_Line\event'
output_dir = 'mslp+tcrw+quivers_output'
os.makedirs(output_dir, exist_ok=True)

# Load datasets
ds_mslp = xr.open_dataset(os.path.join(input_dir, 'mslp.nc'))
ds_tcrw = xr.open_dataset(os.path.join(input_dir, 'tcrw.nc'))
ds_u10 = xr.open_dataset(os.path.join(input_dir, 'u10.nc'))
ds_v10 = xr.open_dataset(os.path.join(input_dir, 'v10.nc'))

mslp = ds_mslp['msl']
tcrw = ds_tcrw['tcrw']
u10 = ds_u10['u10']
v10 = ds_v10['v10']

extent = [100, 150, 0, 60]

def plot_mslp_with_quivers(mslp, u10, v10, time, index, extent, output_dir):
    mslp_t = mslp.sel(valid_time=time).squeeze().values / 100.0  # Pa to hPa
    u10_t = u10.sel(valid_time=time).values
    v10_t = v10.sel(valid_time=time).values

    lon = mslp.longitude.values
    lat = mslp.latitude.values

    lon2d, lat2d = np.meshgrid(lon, lat)

    fig, ax = plt.subplots(figsize=(12, 6), subplot_kw={'projection': ccrs.PlateCarree()})
    ax.set_extent(extent, crs=ccrs.PlateCarree())
    ax.coastlines(resolution='10m')
    ax.add_feature(cfeature.BORDERS.with_scale('10m'), linewidth=0.5)
    ax.add_feature(cfeature.LAND.with_scale('10m'), facecolor='lightgray')
    ax.add_feature(cfeature.RIVERS.with_scale('10m'))

    pcm = ax.pcolormesh(
        lon2d, lat2d, mslp_t,
        cmap='coolwarm_r',
        vmin=980,
        vmax=1040,
        shading='auto',
        transform=ccrs.PlateCarree()
    )

    cbar = plt.colorbar(pcm, ax=ax, orientation='horizontal', pad=0.05)
    cbar.set_label('Mean Sea Level Pressure (hPa)')

    skip = 5
    ax.quiver(
        u10.longitude.values[::skip], u10.latitude.values[::skip],
        u10_t[::skip, ::skip], v10_t[::skip, ::skip],
        scale=700, width=0.002, color='black',
        transform=ccrs.PlateCarree()
    )

    timestamp = np.datetime_as_string(time, unit='h')
    ax.set_title(f"MSLP with Wind Quivers - {timestamp} UTC")

    plt.savefig(os.path.join(output_dir, f"mslp_quiver_{index:03}.png"), dpi=300, bbox_inches='tight')
    plt.close()

def plot_tcrw_with_quivers(tcrw, u10, v10, time, index, extent, output_dir):
    tcrw_t = tcrw.sel(valid_time=time).values.squeeze()
    u10_t = u10.sel(valid_time=time).values
    v10_t = v10.sel(valid_time=time).values

    lon = tcrw.longitude.values
    lat = tcrw.latitude.values

    lon2d, lat2d = np.meshgrid(lon, lat)

    fig, ax = plt.subplots(figsize=(12, 6), subplot_kw={'projection': ccrs.PlateCarree()})
    ax.set_extent(extent, crs=ccrs.PlateCarree())
    ax.coastlines(resolution='10m')
    ax.add_feature(cfeature.BORDERS.with_scale('10m'), linewidth=0.5)
    ax.add_feature(cfeature.LAND.with_scale('10m'), facecolor='lightgray')
    ax.add_feature(cfeature.RIVERS.with_scale('10m'))

    pcm = ax.pcolormesh(
        lon2d, lat2d, tcrw_t,
        cmap='Reds',
        vmin=0,
        vmax=2,
        shading='auto',
        transform=ccrs.PlateCarree()
    )
    cbar = plt.colorbar(pcm, ax=ax, orientation='horizontal', pad=0.05)
    cbar.set_label('Total Column Rain Water (kg/m² ≈ mm)')

    skip = 5
    ax.quiver(
        u10.longitude.values[::skip], u10.latitude.values[::skip],
        u10_t[::skip, ::skip], v10_t[::skip, ::skip],
        scale=700, width=0.002, color='black',
        transform=ccrs.PlateCarree()
    )

    timestamp = np.datetime_as_string(time, unit='h')
    ax.set_title(f"TCRW with Wind Quivers - {timestamp} UTC")

    filename = os.path.join(output_dir, f"tcrw_quiver_{index:03}.png")
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()

if __name__ == '__main__':
    # Find common times across all datasets
    valid_times = np.intersect1d(
        np.intersect1d(mslp.valid_time.values, tcrw.valid_time.values),
        np.intersect1d(u10.valid_time.values, v10.valid_time.values)
    )

    for i, time in enumerate(valid_times):
        print(f"Plotting time index {i}: {np.datetime_as_string(time, unit='h')}")

        plot_mslp_with_quivers(mslp, u10, v10, time, i, extent, output_dir)
        plot_tcrw_with_quivers(tcrw, u10, v10, time, i, extent, output_dir)

    ds_mslp.close()
    ds_tcrw.close()
    ds_u10.close()
    ds_v10.close()
