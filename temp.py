import os
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# ===== USER SETTINGS =====
input_file = r'mean_wind_30years_test_output\3yr_NDJFM_mean_wind.nc'  # <- update path if needed
output_dir = 'temp_output'
os.makedirs(output_dir, exist_ok=True)

extent = [100, 150, 0, 60]  # [lon_min, lon_max, lat_min, lat_max]
lon_min, lon_max, lat_min, lat_max = extent

# ===== PLOTTING FUNCTIONS =====
def plot_streamlines(u, v, lon, lat, title, save_path):
    wind_speed = np.sqrt(u**2 + v**2)
    lon_grid, lat_grid = np.meshgrid(lon, lat)

    cmap_colors = ["#dedede", "#95cee2", "#c1e4c4", "#bbdc71", "#8cce39", "#0eab42",
                   "#fbee6b", "#eecf52", "#eea056", "#db4218"]
    cmap_levels = [0, 1, 2, 4, 6, 9, 11, 14, 17, 21, 24]
    cmap = mcolors.ListedColormap(cmap_colors)
    norm = mcolors.BoundaryNorm(cmap_levels, cmap.N)

    fig, ax = plt.subplots(figsize=(10, 6), subplot_kw={"projection": ccrs.PlateCarree()})
    ax.set_extent(extent, crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    ax.add_feature(cfeature.LAND, edgecolor='black', facecolor='lightgray')

    gl = ax.gridlines(draw_labels=True, linestyle="--", linewidth=0.5)
    gl.top_labels = False
    gl.right_labels = False

    contour = ax.pcolormesh(lon_grid, lat_grid, wind_speed, cmap=cmap, norm=norm, transform=ccrs.PlateCarree())
    ax.streamplot(lon_grid, lat_grid, u, v, density=2.0, color='black', transform=ccrs.PlateCarree())

    fig.subplots_adjust(top=0.85)
    plt.figtext(0.5, 0.97, title, ha='center', va='bottom', fontsize=14, fontweight='bold')

    cbar = plt.colorbar(contour, ax=ax, orientation="vertical", pad=0.02, aspect=25)
    cbar.set_label("Wind Speed (m/s)")
    cbar.set_ticks(cmap_levels)

    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"Streamline plot saved: {save_path}")
    plt.close(fig)

def plot_quiver(u, v, lon, lat, title, save_path):
    fig, ax = plt.subplots(figsize=(10, 6), subplot_kw={"projection": ccrs.PlateCarree()})
    ax.set_extent(extent, crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    ax.add_feature(cfeature.LAND, edgecolor='black', facecolor='lightgray')

    gl = ax.gridlines(draw_labels=True, linestyle="--", linewidth=0.5)
    gl.top_labels = False
    gl.right_labels = False

    lon_grid, lat_grid = np.meshgrid(lon, lat)
    skip = (slice(None, None, 5), slice(None, None, 5))  # skip every 5th arrow to declutter
    ax.quiver(
        lon_grid[skip], lat_grid[skip],
        u[skip], v[skip],
        transform=ccrs.PlateCarree(),
        scale=200, width=0.005, color='blue'
    )

    fig.subplots_adjust(top=0.85)
    plt.figtext(0.5, 0.97, title, ha='center', va='bottom', fontsize=14, fontweight='bold')

    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"Quiver plot saved: {save_path}")
    plt.close(fig)

# ===== MAIN SCRIPT =====
# Load the NetCDF file
ds = xr.open_dataset(input_file)

# Extract variables
u_mean = ds['u10_mean'].values
v_mean = ds['v10_mean'].values
lon = ds['longitude'].values
lat = ds['latitude'].values

# Generate plots
streamline_path = os.path.join(output_dir, 'mean_streamline.png')
quiver_path = os.path.join(output_dir, 'mean_quiver.png')

plot_streamlines(u_mean, v_mean, lon, lat, "Mean 10-m Streamlines", streamline_path)
plot_quiver(u_mean, v_mean, lon, lat, "Mean 10-m Wind Vectors", quiver_path)

print("All plots generated successfully!")
