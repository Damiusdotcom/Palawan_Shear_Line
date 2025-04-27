def save_mean_to_netcdf(u_mean, v_mean, lon, lat, save_path):
    """Save the mean u and v wind components into a NetCDF file with metadata."""
    ds = xr.Dataset(
        {
            "u10_mean": (["latitude", "longitude"], u_mean.values, {
                "long_name": "30-year mean 10-m zonal wind component",
                "units": "m s-1"
            }),
            "v10_mean": (["latitude", "longitude"], v_mean.values, {
                "long_name": "30-year mean 10-m meridional wind component",
                "units": "m s-1"
            }),
        },
        coords={
            "longitude": ("longitude", lon, {"units": "degrees_east"}),
            "latitude": ("latitude", lat, {"units": "degrees_north"}),
        },
        attrs={
            "description": f"{num_years}-year NDJFM mean 10-m wind components",
            "created_by": "Your Name or Script",
            "note": "ND from previous year, JFM from current year.",
            "conventions": "CF-1.6",
        }
    )
    ds.to_netcdf(save_path)
    print(f"Mean wind components saved to NetCDF: {save_path}")


# After calculating the 30-year (or however many years) mean u and v components:
num_years = len(years)

# Set the dynamic filename for the output NetCDF file
final_netcdf_path = os.path.join(output_dir, f'{num_years}yr_NDJFM_mean_wind.nc')

# Call the function to save the data as NetCDF
save_mean_to_netcdf(final_u, final_v, lon, lat, final_netcdf_path)
