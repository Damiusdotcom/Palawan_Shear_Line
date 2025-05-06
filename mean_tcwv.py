import xarray as xr
import numpy as np
import os
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from glob import glob
from dask.diagnostics import ProgressBar

def compute_tcrw_mean_dask(files):
    total_sum = None
    total_count = 0

    for file in files:
        ds = xr.open_dataset(file, chunks={'valid_time': 10})  # Dask-friendly chunking
        tcrw = ds['tcrw']  # Extract total column rain water

        sum_this_file = tcrw.sum(dim='valid_time', skipna=True)
        count_this_file = tcrw['valid_time'].size

        if total_sum is None:
            total_sum = sum_this_file
        else:
            total_sum += sum_this_file

        total_count += count_this_file

    mean_tcrw = total_sum / total_count

    with ProgressBar():
        return mean_tcrw.compute()

def plot_tcrw(tcrw_mean, output_path):
    fig = plt.figure(figsize=(12, 6))
    ax = plt.axes(projection=ccrs.PlateCarree())

    # Set domain extent for Southeast Asia
    extent = [100, 150, 0, 60]
    ax.set_extent(extent, crs=ccrs.PlateCarree())

    ax.coastlines()
    ax.add_feature(cfeature.BORDERS, linewidth=0.5)
    ax.add_feature(cfeature.LAND, facecolor='lightgray')
    ax.add_feature(cfeature.LAKES, edgecolor='black', facecolor='none')
    ax.add_feature(cfeature.RIVERS)

    contour = ax.contourf(tcrw_mean.longitude, tcrw_mean.latitude, tcrw_mean,
                          levels=np.arange(0, 150, 5), cmap='Blues', extend='both')

    cbar = plt.colorbar(contour, ax=ax, orientation='horizontal', pad=0.05)
    cbar.set_label('Total Column Rain Water (kg/m²)')

    ax.set_title('Mean Total Column Rain Water (1991–2020)', fontsize=14)
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

def save_tcrw_mean_netcdf(tcrw_mean, output_path):
    tcrw_mean.name = 'tcrw'
    tcrw_mean.attrs['units'] = 'kg/m²'
    tcrw_mean.attrs['long_name'] = 'Mean Total Column Rain Water'
    tcrw_mean.to_netcdf(output_path)

if __name__ == "__main__":
    input_folder = r'Y:\wfs_shared\Personnel Files\DJEV\Palawan_Shear_Line\TCWV'  # folder with all yearly NetCDF files
    output_folder = "mean_tcwv_output"
    os.makedirs(output_folder, exist_ok=True)

    output_plot = os.path.join(output_folder, "mean_tcrw_plot.png")
    output_netcdf = os.path.join(output_folder, "mean_tcrw_mean.nc")

    nc_files = sorted(glob(os.path.join(input_folder, "*.nc")))

    if not nc_files:
        raise FileNotFoundError("No NetCDF files found in the specified folder.")

    print(f"Processing {len(nc_files)} files...")

    tcrw_mean = compute_tcrw_mean_dask(nc_files)
    save_tcrw_mean_netcdf(tcrw_mean, output_netcdf)
    plot_tcrw(tcrw_mean, output_plot)

    print("✅ Mean TCWV NetCDF and plot saved successfully.")
