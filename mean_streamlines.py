import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import os
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D

# File paths
u_file_2019 = 'Y:/wfs_shared/Personnel Files/DJEV/Palawan_Shear_Line/mean_streamlines/u_2019.nc'
v_file_2019 = 'Y:/wfs_shared/Personnel Files/DJEV/Palawan_Shear_Line/mean_streamlines/v_2019.nc'
u_file_2020 = 'Y:/wfs_shared/Personnel Files/DJEV/Palawan_Shear_Line/mean_streamlines/u_2020.nc'
v_file_2020 = 'Y:/wfs_shared/Personnel Files/DJEV/Palawan_Shear_Line/mean_streamlines/v_2020.nc'

# Load datasets
ds_u_2019 = xr.open_dataset(u_file_2019)
ds_v_2019 = xr.open_dataset(v_file_2019)
ds_u_2020 = xr.open_dataset(u_file_2020)
ds_v_2020 = xr.open_dataset(v_file_2020)

# Domain
lat_min, lat_max = 0, 60
lon_min, lon_max = 100, 150

# NDJFM months
months_info = {
    "2019-11": (ds_u_2019, ds_v_2019),
    "2019-12": (ds_u_2019, ds_v_2019),
    "2020-01": (ds_u_2020, ds_v_2020),
    "2020-02": (ds_u_2020, ds_v_2020),
    "2020-03": (ds_u_2020, ds_v_2020),
}

output_dir = "mean_streamlines_output/2019_NDJFM"
os.makedirs(output_dir, exist_ok=True)

# Custom colormap for streamlines
cmap_colors = ["#dedede", "#95cee2", "#c1e4c4", "#bbdc71", "#8cce39", "#0eab42",
               "#fbee6b", "#eecf52", "#eea056", "#db4218"]
cmap_levels = [0, 1, 2, 4, 6, 9, 11, 14, 17, 21, 24]
cmap = mcolors.ListedColormap(cmap_colors)
norm = mcolors.BoundaryNorm(cmap_levels, cmap.N)

monthly_u = []
monthly_v = []

for month, (ds_u, ds_v) in months_info.items():
    u_month = ds_u.sel(valid_time=ds_u.valid_time.dt.month == int(month[-2:]))
    v_month = ds_v.sel(valid_time=ds_v.valid_time.dt.month == int(month[-2:]))

    u_mean = u_month.u10.mean(dim='valid_time')
    v_mean = v_month.v10.mean(dim='valid_time')

    lat = u_mean.latitude.sel(latitude=slice(lat_max, lat_min)).values
    lon = u_mean.longitude.sel(longitude=slice(lon_min, lon_max)).values

    u_mean = u_mean.sel(latitude=slice(lat_max, lat_min), longitude=slice(lon_min, lon_max)).values
    v_mean = v_mean.sel(latitude=slice(lat_max, lat_min), longitude=slice(lon_min, lon_max)).values
    wind_speed = np.sqrt(u_mean**2 + v_mean**2)

    lon_grid, lat_grid = np.meshgrid(lon, lat)

    # --- Streamline Plot ---
    fig1, ax1 = plt.subplots(figsize=(10, 6), subplot_kw={"projection": ccrs.PlateCarree()})
    ax1.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    ax1.add_feature(cfeature.COASTLINE)
    ax1.add_feature(cfeature.BORDERS, linestyle=':')
    ax1.add_feature(cfeature.LAND, edgecolor='black', facecolor='lightgray')
    gl = ax1.gridlines(draw_labels=True, linestyle="--", linewidth=0.5)
    gl.top_labels = False
    gl.right_labels = False

    contour = ax1.pcolormesh(lon, lat, wind_speed, cmap=cmap, norm=norm, transform=ccrs.PlateCarree())
    ax1.streamplot(lon_grid, lat_grid, u_mean, v_mean, density=2.0, color='black', transform=ccrs.PlateCarree())

    title = f"10 m Streamlines - {month}"
    fig1.subplots_adjust(top=0.85)
    plt.figtext(0.5, 0.97, title, ha='center', va='bottom', fontsize=14, fontweight='bold', color='black')

    cbar = plt.colorbar(contour, ax=ax1, orientation="vertical", pad=0.02, aspect=25)
    cbar.set_label("Wind Speed (m/s)")
    cbar.set_ticks(cmap_levels)

    save_path = f"{output_dir}/streamline_10m_{month}.png"
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"Streamline plot saved: {save_path}")
    plt.close(fig1)

    # --- Quiver Plot ---
    fig2, ax2 = plt.subplots(figsize=(10, 6), subplot_kw={"projection": ccrs.PlateCarree()})
    ax2.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    ax2.add_feature(cfeature.COASTLINE)
    ax2.add_feature(cfeature.BORDERS, linestyle=':')
    ax2.add_feature(cfeature.LAND, edgecolor='black', facecolor='lightgray')
    gl = ax2.gridlines(draw_labels=True, linestyle="--", linewidth=0.5)
    gl.top_labels = False
    gl.right_labels = False

    # Shorten the quivers and make them thicker
    skip = (slice(None, None, 5), slice(None, None, 5))
    q = ax2.quiver(lon_grid[skip], lat_grid[skip], u_mean[skip], v_mean[skip], 
                   transform=ccrs.PlateCarree(), scale=200, width=0.005, color='blue')

    # Move arrow to higher and more left position
    ref_lon = lon.min() + 2  # further left (increase to move right)
    ref_lat = lat.max() - 0.5  # higher (decrease to move even higher)

    # 2. Draw the reference arrow (pointing east)
    ax2.quiver(ref_lon, ref_lat, 10, 0,
           transform=ccrs.PlateCarree(),
           scale=200, width=0.005, color='blue')

    # 3. Add a label right next to the arrow
    ax2.text(ref_lon + 1.5, ref_lat, '10 m/s',
         transform=ccrs.PlateCarree(),
         fontsize=10, va='center')

    # Adjust layout to ensure the key isn’t cut off
    plt.tight_layout(rect=[0, 0.05, 1, 1])

    title = f"10 m Wind Vectors (Quiver) - {month}"
    fig2.subplots_adjust(top=0.85)
    plt.figtext(0.5, 0.97, title, ha='center', va='bottom', fontsize=14, fontweight='bold', color='black')

    save_path = f"{output_dir}/quiver_10m_{month}.png"
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"Quiver plot saved: {save_path}")
    plt.close(fig2)

    monthly_u.append(u_mean)
    monthly_v.append(v_mean)

# --- NDJFM MEAN ---
u_ndjfm = np.mean(monthly_u, axis=0)
v_ndjfm = np.mean(monthly_v, axis=0)
wind_ndjfm = np.sqrt(u_ndjfm**2 + v_ndjfm**2)

# Streamline NDJFM
fig1, ax1 = plt.subplots(figsize=(10, 6), subplot_kw={"projection": ccrs.PlateCarree()})
ax1.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
ax1.add_feature(cfeature.COASTLINE)
ax1.add_feature(cfeature.BORDERS, linestyle=':')
ax1.add_feature(cfeature.LAND, edgecolor='black', facecolor='lightgray')
gl = ax1.gridlines(draw_labels=True, linestyle="--", linewidth=0.5)
gl.top_labels = False
gl.right_labels = False

contour = ax1.pcolormesh(lon, lat, wind_ndjfm, cmap=cmap, norm=norm, transform=ccrs.PlateCarree())
ax1.streamplot(lon_grid, lat_grid, u_ndjfm, v_ndjfm, density=2.0, color='black', transform=ccrs.PlateCarree())

title = "10 m Streamlines - 2019 NDJFM Mean"
fig1.subplots_adjust(top=0.85)
plt.figtext(0.5, 0.97, title, ha='center', va='bottom', fontsize=14, fontweight='bold', color='black')

cbar = plt.colorbar(contour, ax=ax1, orientation="vertical", pad=0.02, aspect=25)
cbar.set_label("Wind Speed (m/s)")
cbar.set_ticks(cmap_levels)

save_path = f"{output_dir}/streamline_10m_2019_NDJFM.png"
plt.savefig(save_path, dpi=300, bbox_inches='tight')
print(f"Streamline NDJFM plot saved: {save_path}")
plt.close(fig1)

# Quiver NDJFM
fig2, ax2 = plt.subplots(figsize=(10, 6), subplot_kw={"projection": ccrs.PlateCarree()})
ax2.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
ax2.add_feature(cfeature.COASTLINE)
ax2.add_feature(cfeature.BORDERS, linestyle=':')
ax2.add_feature(cfeature.LAND, edgecolor='black', facecolor='lightgray')
gl = ax2.gridlines(draw_labels=True, linestyle="--", linewidth=0.5)
gl.top_labels = False
gl.right_labels = False

# Shorten the quivers and make them thicker
skip = (slice(None, None, 5), slice(None, None, 5))
q = ax2.quiver(lon_grid[skip], lat_grid[skip], u_ndjfm[skip], v_ndjfm[skip],
               transform=ccrs.PlateCarree(), scale=200, width=0.005, color='blue')

# Reference arrow: place it inside the axes (e.g., bottom right corner)
ax2.quiverkey(q, X=0.85, Y=0.05, U=10,
              label='10 m/s', labelpos='E', coordinates='axes')

# Adjust layout to ensure the key isn’t cut off
plt.tight_layout(rect=[0, 0.05, 1, 1])

title = "10 m Wind Vectors (Quiver) - 2019 NDJFM Mean"
fig2.subplots_adjust(top=0.85)
plt.figtext(0.5, 0.97, title, ha='center', va='bottom', fontsize=14, fontweight='bold', color='black')

save_path = f"{output_dir}/quiver_10m_2019_NDJFM.png"
plt.savefig(save_path, dpi=300, bbox_inches='tight')
print(f"Quiver NDJFM plot saved: {save_path}")
plt.close(fig2)
