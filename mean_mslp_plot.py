import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import os

def plot_mslp_with_quivers(mslp_file, wind_file, extent, output_path, title):
  
    lon_min, lon_max, lat_min, lat_max = extent

    # Load and clip MSLP
    ds_mslp = xr.open_dataset(mslp_file)
    ds_mslp = ds_mslp.sel(latitude=slice(lat_max, lat_min), longitude=slice(lon_min, lon_max))
    mslp = ds_mslp['msl']

    # Load and clip wind data
    ds_wind = xr.open_dataset(wind_file)
    ds_wind = ds_wind.sel(latitude=slice(lat_max, lat_min), longitude=slice(lon_min, lon_max))
    u = ds_wind['u10_mean']
    v = ds_wind['v10_mean']

    # Create plot
    fig, ax = plt.subplots(figsize=(10, 6), subplot_kw={"projection": ccrs.PlateCarree()})
    ax.set_extent(extent)
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    ax.add_feature(cfeature.LAND, edgecolor='black', facecolor='lightgray')
    ax.add_feature(cfeature.LAKES, edgecolor='black', facecolor='none')
    ax.add_feature(cfeature.RIVERS)

    gl = ax.gridlines(draw_labels=True, linestyle="--", linewidth=0.5)
    gl.top_labels = False
    gl.right_labels = False

    msl_hpa = mslp / 100  # Convert Pa to hPa

    # TCWV contour plot
    contour = ax.contourf(ds_mslp.longitude, ds_mslp.latitude, msl_hpa,
                          levels=np.arange(980, 1040, 2), cmap='coolwarm_r', extend='both')

    # Quiver overlay for wind vectors
    lon_grid, lat_grid = np.meshgrid(ds_wind.longitude, ds_wind.latitude)
    skip = (slice(None, None, 5), slice(None, None, 5))  # Downsample the wind grid for clarity
    ax.quiver(
        lon_grid[skip], lat_grid[skip],
        u.values[skip], v.values[skip],
        transform=ccrs.PlateCarree(), scale=200, width=0.005, color='red'
    )

    # Colorbar
    cbar = plt.colorbar(contour, ax=ax, orientation='horizontal', pad=0.05)
    cbar.set_label('MSLP (hPa)')

    # Title and save
    ax.set_title(title, fontsize=14)
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"MSLP + Quiver plot saved: {output_path}")
    plt.close()


# ======= CONFIGURATION =======
mslp_nc_file = "mean_mslp_output/ndjfm_mslp_mean.nc"
wind_nc_file = "mean_wind_30years_output/30yr_NDJFM_mean_wind.nc"

extent = [100, 150, 0, 60]  # Define the area of interest (Southeast Asia)
output_dir = "mean_mslp_plot_outputs"
os.makedirs(output_dir, exist_ok=True)

output_file = os.path.join(output_dir, "mean_mslp_with_wind_vectors.png")
plot_title = "Mean MSLP with 10m Wind Vectors (1991–2020)"

# Call the function
plot_mslp_with_quivers(mslp_nc_file, wind_nc_file, extent, output_file, plot_title)
