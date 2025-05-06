import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import os
import pandas as pd

def plot_anomalies_with_quivers(variable_file, mean_file, u10_file, v10_file, wind_mean_file,
                                 extent, output_dir, title_prefix, variable_name, date_range, units=""):
    """
    Plots anomalies of a variable (MSLP or TCWV) with 10-m wind anomaly quiver overlay.

    Parameters:
        variable_file (str): Path to NetCDF file containing daily data for the variable (MSLP or TCWV).
        mean_file (str): Path to NetCDF file containing 30-year mean data for the variable (MSLP or TCWV).
        u10_file (str): Path to NetCDF file containing 10-day u10 data.
        v10_file (str): Path to NetCDF file containing 10-day v10 data.
        wind_mean_file (str): Path to NetCDF file containing 'u10_mean' and 'v10_mean' for 30-year average.
        extent (list): [lon_min, lon_max, lat_min, lat_max].
        output_dir (str): Directory to save the output plots.
        title_prefix (str): Title prefix for the plot.
        variable_name (str): Name of the variable (either 'mslp' or 'tcrw').
        date_range (list): List of dates corresponding to the 10-day period.
        units (str): Units of the variable (default is an empty string).
    """
    lon_min, lon_max, lat_min, lat_max = extent

    # Load and clip the 30-year mean variable
    ds_mean = xr.open_dataset(mean_file).sel(latitude=slice(lat_max, lat_min), longitude=slice(lon_min, lon_max))
    mean_var = ds_mean[variable_name]

    # Load and clip the 10-day period variable
    ds_variable = xr.open_dataset(variable_file).sel(latitude=slice(lat_max, lat_min), longitude=slice(lon_min, lon_max))
    var = ds_variable[variable_name]

    # Load and clip 10-day u and v wind data
    ds_u = xr.open_dataset(u10_file).sel(latitude=slice(lat_max, lat_min), longitude=slice(lon_min, lon_max))
    ds_v = xr.open_dataset(v10_file).sel(latitude=slice(lat_max, lat_min), longitude=slice(lon_min, lon_max))
    u10 = ds_u['u10']
    v10 = ds_v['v10']

    # Load and clip 30-year mean wind data
    ds_wind_mean = xr.open_dataset(wind_mean_file).sel(latitude=slice(lat_max, lat_min), longitude=slice(lon_min, lon_max))
    u10_mean = ds_wind_mean['u10_mean']
    v10_mean = ds_wind_mean['v10_mean']

    # Compute anomalies
    var_anomaly = var - mean_var
    u_anomaly = u10 - u10_mean
    v_anomaly = v10 - v10_mean

    # Define constant color levels for the contour plot
    var_min = np.min(var_anomaly.values)
    var_max = np.max(var_anomaly.values)
    color_levels = np.linspace(var_min, var_max, 21)

    # Plot for each day
    for i in range(var_anomaly.shape[0]):
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

        # Plot contour for the anomaly of the variable with constant color levels
        contour = ax.contourf(
            ds_variable.longitude, ds_variable.latitude, var_anomaly[i].values,
            levels=color_levels,
            cmap='coolwarm_r', extend='both'
        )

        # Quiver overlay for wind anomalies
        lon_grid, lat_grid = np.meshgrid(ds_u.longitude, ds_u.latitude)
        skip = (slice(None, None, 5), slice(None, None, 5))  # To reduce the density of the arrows
        ax.quiver(
            lon_grid[skip], lat_grid[skip],
            u_anomaly[i].values[skip], v_anomaly[i].values[skip],
            transform=ccrs.PlateCarree(), scale=200, width=0.005, color='black'
        )

        # Add colorbar and title
        cbar = plt.colorbar(contour, ax=ax, orientation='horizontal', pad=0.05)
        cbar.set_label(f'{variable_name.upper()} Anomaly ({units})')

        # Calculate day label and actual date
        day_offset = i - 5  # Day -5 to Day +5 relative to Day 0
        actual_date = pd.to_datetime(date_range[5]) + pd.Timedelta(days=day_offset)
        title = f"{title_prefix} - Day {day_offset} ({actual_date.strftime('%Y-%m-%d')})"
        ax.set_title(title, fontsize=14)

        # Save the plot with only the date in the file name (no day offset)
        output_file = os.path.join(output_dir, f"{variable_name}_{actual_date.strftime('%Y-%m-%d')}_anomaly.png")
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Plot saved: {output_file} for {variable_name} Day {day_offset}")

    print(f"✅ Finished processing {variable_name} anomalies from: {variable_file}")

# ======= CONFIGURATION =======

# 30-year mean files
mslp_mean_file = "mean_mslp_output/ndjfm_mslp_mean.nc"
tcrw_mean_file = "mean_tcwv_output/mean_tcrw_mean.nc"
wind_mean_file = "mean_wind_30years_output/30yr_NDJFM_mean_wind.nc"

# 10-day period files (separate u and v)
mslp_10day_file = r'Y:\wfs_shared\Personnel Files\DJEV\Palawan_Shear_Line\event\mslp.nc'
tcrw_10day_file = r'Y:\wfs_shared\Personnel Files\DJEV\Palawan_Shear_Line\event\tcrw.nc'
u10_10day_file = r'Y:\wfs_shared\Personnel Files\DJEV\Palawan_Shear_Line\event\u10.nc'
v10_10day_file = r'Y:\wfs_shared\Personnel Files\DJEV\Palawan_Shear_Line\event\v10.nc'

# Area of interest
extent = [100, 150, 0, 60]

# Output directory
output_dir = "mslp+tcwv_anomaly_outputs"
os.makedirs(output_dir, exist_ok=True)

# Titles
mslp_title = "MSLP Anomaly with Wind Anomaly Vectors"
tcrw_title = "TCWV Anomaly with Wind Anomaly Vectors"

# Date range for the 10-day period (should be in the format 'YYYY-MM-DD')
date_range = ['2025-02-04', '2025-02-05', '2025-02-06', '2025-02-07', '2025-02-08', 
              '2025-02-09', '2025-02-10', '2025-02-11', '2025-02-12', '2025-02-13', '2025-02-14']

# Generate anomaly plots
plot_anomalies_with_quivers(mslp_10day_file, mslp_mean_file, u10_10day_file, v10_10day_file, wind_mean_file,
                             extent, output_dir, mslp_title, 'msl', date_range, units="hPa")

plot_anomalies_with_quivers(tcrw_10day_file, tcrw_mean_file, u10_10day_file, v10_10day_file, wind_mean_file,
                             extent, output_dir, tcrw_title, 'tcrw', date_range, units="mm")

print("✅ All anomaly plots generated successfully.")

