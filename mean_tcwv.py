import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from glob import glob
from dask.diagnostics import ProgressBar

def compute_tcwv_mean_daily_dask(files):
    total_sum = None
    total_count = 0

    for file in files:
        ds = xr.open_dataset(file)
        ds = ds.chunk({'valid_time': 240})  # hourly data, daily means

        tcwv = ds['tcwv']

        # Resample to daily mean first
        tcwv_daily = tcwv.resample(valid_time='1D').mean()

        # Decode time for filtering
        time = tcwv_daily.valid_time

        # Filter for ND (1991) and JFM (2021)
        time_filtered = time[
            ((time.dt.year == 1991) & (time.dt.month.isin([11, 12]))) |
            ((time.dt.year == 2021) & (time.dt.month.isin([1, 2, 3])))
        ]

        if time_filtered.size == 0:
            ds.close()
            continue  # skip if no dates

        # Select data
        tcwv_filtered = tcwv_daily.sel(valid_time=time_filtered)

        sum_this_file = tcwv_filtered.sum(dim='valid_time', skipna=True)
        count_this_file = tcwv_filtered.sizes['valid_time']

        if total_sum is None:
            total_sum = sum_this_file
        else:
            total_sum += sum_this_file

        total_count += count_this_file

        ds.close()

    if total_count == 0:
        raise ValueError("No matching NDJFM daily dates found in the data.")

    mean_tcwv = total_sum / total_count

    with ProgressBar():
        return mean_tcwv.compute()

def plot_tcwv(tcwv_mean, output_path):
    fig = plt.figure(figsize=(12, 6))
    ax = plt.axes(projection=ccrs.PlateCarree())

    extent = [100, 150, 0, 60]
    ax.set_extent(extent, crs=ccrs.PlateCarree())

    ax.coastlines()
    ax.add_feature(cfeature.BORDERS, linewidth=0.5)
    ax.add_feature(cfeature.LAND, facecolor='lightgray')
    ax.add_feature(cfeature.LAKES, edgecolor='black', facecolor='none')
    ax.add_feature(cfeature.RIVERS)

    contour = ax.contourf(tcwv_mean.longitude, tcwv_mean.latitude, tcwv_mean,
                          levels=np.arange(0, 70, 2), cmap='YlGnBu', extend='both')

    cbar = plt.colorbar(contour, ax=ax, orientation='horizontal', pad=0.05)
    cbar.set_label('Total Column Water Vapor (kg/m²)')

    ax.set_title('30-Year Mean Daily Total Column Water Vapor (NDJFM)', fontsize=14)
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

def save_tcwv_mean_netcdf(tcwv_mean, output_path):
    tcwv_mean.name = 'tcwv'
    tcwv_mean.attrs['units'] = 'kg/m²'
    tcwv_mean.attrs['long_name'] = 'Mean Daily Total Column Water Vapor'
    tcwv_mean.to_netcdf(output_path)

if __name__ == "__main__":
    input_folder = r'Y:\wfs_shared\Personnel Files\DJEV\Palawan_Shear_Line\TCWV'
    output_folder = "mean_tcwv_output"
    os.makedirs(output_folder, exist_ok=True)

    output_plot = os.path.join(output_folder, "mean_tcwv_plot.png")
    output_netcdf = os.path.join(output_folder, "mean_tcwv_mean.nc")

    nc_files = sorted(glob(os.path.join(input_folder, "*.nc")))

    if not nc_files:
        raise FileNotFoundError("No NetCDF files found in the specified folder.")

    print(f"Processing {len(nc_files)} files...")

    tcwv_mean = compute_tcwv_mean_daily_dask(nc_files)
    save_tcwv_mean_netcdf(tcwv_mean, output_netcdf)
    plot_tcwv(tcwv_mean, output_plot)

    print("✅ Mean daily TCWV NetCDF and plot saved successfully.")
