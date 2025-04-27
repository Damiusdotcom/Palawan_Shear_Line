import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# ========== USER SETTINGS ==========
input_dir = r'Y:\wfs_shared\Personnel Files\DJEV\Palawan_Shear_Line\mean_streamlines'  # Update to your directory
output_dir = 'mean_wind_30years_output'
os.makedirs(output_dir, exist_ok=True)

start_year = 2012
end_year = 2020  # Last year to start NDJFM (will fetch until Mar 2021)
years = np.arange(start_year, end_year + 1)

months = [11, 12, 1, 2, 3]  # Nov to Mar
month_names = ['Nov', 'Dec', 'Jan', 'Feb', 'Mar']

# ========== LOAD ALL DATA ==========
u_files = {}
v_files = {}
for year in range(start_year, end_year + 2):  # need +1 year for JFM of next year
    u_files[year] = xr.open_dataset(f'{input_dir}/u_{year}.nc')
    v_files[year] = xr.open_dataset(f'{input_dir}/v_{year}.nc')

# ========== SET UP DOMAIN ==========
extent = [100, 150, 0, 60]  # [lon_min, lon_max, lat_min, lat_max]
lon_min, lon_max, lat_min, lat_max = extent

lon = u_files[start_year]['longitude'].values
lat = u_files[start_year]['latitude'].values

# ========== STORAGE FOR FINAL MEAN ==========
all_ndjfm_u = []
all_ndjfm_v = []

# ========== HELPER FUNCTIONS ==========
def plot_streamlines(u, v, title, save_path):
    wind_speed = np.sqrt(u.values**2 + v.values**2)
    lon_grid, lat_grid = np.meshgrid(u.coords['longitude'].values, u.coords['latitude'].values)

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
    ax.streamplot(lon_grid, lat_grid, u.values, v.values, density=2.0, color='black', transform=ccrs.PlateCarree())

    fig.subplots_adjust(top=0.85)
    plt.figtext(0.5, 0.97, title, ha='center', va='bottom', fontsize=14, fontweight='bold', color='black')

    cbar = plt.colorbar(contour, ax=ax, orientation="vertical", pad=0.02, aspect=25)
    cbar.set_label("Wind Speed (m/s)")
    cbar.set_ticks(cmap_levels)

    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"Streamline plot saved: {save_path}")
    plt.close(fig)

def plot_quiver(u, v, lon, lat, title, save_path):
    fig, ax = plt.subplots(figsize=(10, 6), subplot_kw={"projection": ccrs.PlateCarree()})
    ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    ax.add_feature(cfeature.LAND, edgecolor='black', facecolor='lightgray')

    gl = ax.gridlines(draw_labels=True, linestyle="--", linewidth=0.5)
    gl.top_labels = False
    gl.right_labels = False

    lon_grid, lat_grid = np.meshgrid(lon, lat)
    skip = (slice(None, None, 5), slice(None, None, 5))
    ax.quiver(
        lon_grid[skip], lat_grid[skip],
        u.values[skip], v.values[skip],
        transform=ccrs.PlateCarree(),
        scale=200, width=0.005, color='blue'
    )

    plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.05)
    plt.tight_layout()
    plt.figtext(0.5, 0.97, title, ha='center', va='bottom', fontsize=14, fontweight='bold', color='black')

    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"Quiver plot saved: {save_path}")
    plt.close(fig)

# ========== PROCESS EACH YEAR ==========
for year in years:
    print(f'Processing {year} NDJFM...')
    year_folder = os.path.join(output_dir, str(year))
    os.makedirs(year_folder, exist_ok=True)

    monthly_means_u = []
    monthly_means_v = []

    for month in months:
        if month in [11, 12]:
            u_data = u_files[year]
            v_data = v_files[year]
        else:  # Jan-Feb-Mar from next year
            u_data = u_files[year + 1]
            v_data = v_files[year + 1]

        u_month = u_data['u10'].sel(valid_time=u_data['valid_time'].dt.month == month).mean(dim='valid_time')
        v_month = v_data['v10'].sel(valid_time=v_data['valid_time'].dt.month == month).mean(dim='valid_time')

        monthly_means_u.append(u_month)
        monthly_means_v.append(v_month)

        # Adjust year for Jan-Feb-Mar
        if month in [1, 2, 3]:
            plot_year = year + 1
        else:
            plot_year = year

        # Streamline plot for individual month
        streamline_path = os.path.join(year_folder, f'{plot_year}-{month:02d}_streamline.png')
        plot_streamlines(u_month, v_month, f'{plot_year}-{month:02d} Mean 10-m Streamlines', streamline_path)

        # Quiver plot for individual month
        quiver_path = os.path.join(year_folder, f'{plot_year}-{month:02d}_quiver.png')
        plot_quiver(u_month, v_month, lon, lat, f'{plot_year}-{month:02d} Mean 10-m Wind Vectors', quiver_path)


    # Calculate NDJFM mean for the season
    ndjfm_u = sum(monthly_means_u) / len(monthly_means_u)
    ndjfm_v = sum(monthly_means_v) / len(monthly_means_v)

    all_ndjfm_u.append(ndjfm_u)
    all_ndjfm_v.append(ndjfm_v)

    # NDJFM Streamline and Quiver
    ndjfm_streamline_path = os.path.join(year_folder, 'NDJFM_streamline.png')
    plot_streamlines(ndjfm_u, ndjfm_v, f'{year} NDJFM Mean 10-m Streamlines', ndjfm_streamline_path)

    ndjfm_quiver_path = os.path.join(year_folder, 'NDJFM_quiver.png')
    plot_quiver(ndjfm_u, ndjfm_v, lon, lat, f'{year} NDJFM Mean 10-m Wind Vectors', ndjfm_quiver_path)

# ========== FINAL 9-YEAR NDJFM MEAN ==========
final_u = sum(all_ndjfm_u) / len(all_ndjfm_u)
final_v = sum(all_ndjfm_v) / len(all_ndjfm_v)

# After calculating the 30-year (or however many years) mean u and v components:
num_years = len(years)

final_streamline_path = os.path.join(output_dir, f'{num_years}yr_NDJFM_streamline.png')
plot_streamlines(final_u, final_v, f'{num_years} NDJFM Mean 10-m Streamlines', final_streamline_path)

final_quiver_path = os.path.join(output_dir, f'{num_years}yr_NDJFM_quiver.png')
plot_quiver(final_u, final_v, lon, lat, f'{num_years} NDJFM Mean 10-m Wind Vectors', final_quiver_path)

print('All done!')
