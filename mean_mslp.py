import xarray as xr
import numpy as np
import os
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from glob import glob
from dask.diagnostics import ProgressBar

def compute_ndjfm_mean_dask(files):
    total_sum = None
    total_count = 0

    for file in files:
        ds = xr.open_dataset(file)
        ds = ds.chunk({'valid_time': 10})

        # Decode time coordinate to datetime64 for filtering
        time = xr.decode_cf(ds).valid_time

        # Filter for ND 1991 and JFM 2021 months only
        time_filtered = time[
            ((time.dt.year == 1991) & (time.dt.month.isin([11, 12]))) |
            ((time.dt.year == 2021) & (time.dt.month.isin([1, 2, 3])))
        ]

        if time_filtered.size == 0:
            ds.close()
            continue  # skip if no matching dates in this file

        # Select data for filtered times
        msl = ds['msl'].sel(valid_time=time_filtered)

        # Sum over time dimension
        sum_this_file = msl.sum(dim='valid_time', skipna=True)

        # Correct count: get size of filtered time dimension
        count_this_file = msl.sizes['valid_time']

        if total_sum is None:
            total_sum = sum_this_file
        else:
            total_sum += sum_this_file

        total_count += count_this_file

        ds.close()

    if total_count == 0:
        raise ValueError("No matching NDJFM dates found in the data.")

    mean_msl = total_sum / total_count

    with ProgressBar():
        return mean_msl.compute()

def plot_mslp(msl_mean, output_path):
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

    msl_hpa = msl_mean / 100  # Convert Pa to hPa

    contour = ax.contourf(msl_mean.longitude, msl_mean.latitude, msl_hpa,
                          levels=np.arange(980, 1040, 2), cmap='coolwarm_r', extend='both')

    cbar = plt.colorbar(contour, ax=ax, orientation='horizontal', pad=0.05)
    cbar.set_label('Mean Sea Level Pressure (hPa)')

    ax.set_title('NDJFM Mean Sea Level Pressure (1991–2020)', fontsize=14)
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()


def save_ndjfm_mean_netcdf(msl_mean, output_path):
    msl_mean.name = 'msl'
    msl_mean.attrs['units'] = 'Pa'
    msl_mean.attrs['long_name'] = 'NDJFM Mean Sea Level Pressure'
    msl_mean.to_netcdf(output_path)

if __name__ == "__main__":
    input_folder = r'Y:\wfs_shared\Personnel Files\DJEV\Palawan_Shear_Line\MSLP'  # folder with all yearly NetCDF files
    output_folder = "mean_mslp_output"
    os.makedirs(output_folder, exist_ok=True)

    output_plot = os.path.join(output_folder, "ndjfm_mslp_plot.png")
    output_netcdf = os.path.join(output_folder, "ndjfm_mslp_mean.nc")

    nc_files = sorted(glob(os.path.join(input_folder, "*.nc")))

    if not nc_files:
        raise FileNotFoundError("No NetCDF files found in the specified folder.")

    print(f"Processing {len(nc_files)} files...")

    msl_ndjfm_mean = compute_ndjfm_mean_dask(nc_files)
    save_ndjfm_mean_netcdf(msl_ndjfm_mean, output_netcdf)
    plot_mslp(msl_ndjfm_mean, output_plot)

    print("✅ NDJFM MSLP mean NetCDF and plot saved successfully.")
