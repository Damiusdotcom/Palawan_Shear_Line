import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import os

def plot_tcrw_with_quivers(tcrw_file, wind_file, extent, output_path, title):
    """
    Plots Total Column Rain Water (TCWV) with 10-m wind quiver overlay.

    Parameters:
        tcrw_file (str): Path to NetCDF file containing 'tcrw'.
        wind_file (str): Path to NetCDF file containing 'u10_mean' and 'v10_mean'.
        extent (list): [lon_min, lon_max, lat_min, lat_max].
        output_path (str): File path to save the output plot.
        title (str): Title for the plot.
    """
    lon_min, lon_max, lat_min, lat_max = extent

    # Load and clip TCWV
    ds_tcrw = xr.open_dataset(tcrw_file)
    ds_tcrw = ds_tcrw.sel(latitude=slice(lat_max, lat_min), longitude=slice(lon_min, lon_max))
    tcrw = ds_tcrw['tcrw']

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

    # TCWV contour plot
    contour = ax.contourf(ds_tcrw.longitude, ds_tcrw.latitude, tcrw,
                          levels=np.arange(0, 80, 2), cmap='Blues', extend='both')

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
    cbar.set_label('Total Column Rain Water (kg/m²)')

    # Title and save
    ax.set_title(title, fontsize=14)
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"TCWV + Quiver plot saved: {output_path}")
    plt.close()


# ======= CONFIGURATION =======
tcrw_nc_file = "mean_tcwv_output/mean_tcrw_mean.nc"
wind_nc_file = "mean_wind_30years_output/30yr_NDJFM_mean_wind.nc"

extent = [100, 150, 0, 60]  # Define the area of interest (Southeast Asia)
output_dir = "tcwv_plot_outputs"
os.makedirs(output_dir, exist_ok=True)

output_file = os.path.join(output_dir, "mean_tcrw_with_wind_vectors.png")
plot_title = "Mean Total Column Rain Water with 10m Wind Vectors (1991–2020)"

# Call the function
plot_tcrw_with_quivers(tcrw_nc_file, wind_nc_file, extent, output_file, plot_title)
