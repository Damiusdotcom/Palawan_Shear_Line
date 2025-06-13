import xarray as xr

# Change this to your NetCDF file path
file_path = r'Y:\wfs_shared\Personnel Files\DJEV\Palawan_Shear_Line\TCWV\1c811de9cef4708de3a5fbe762d03012.nc'

# Open the NetCDF file
ds = xr.open_dataset(file_path)

# Print global attributes
print("\n📄 Global Attributes:")
for attr_name, attr_value in ds.attrs.items():
    print(f"  {attr_name}: {attr_value}")

# Print variables, their dimensions, and attributes
print("\n📦 Variables:")
for var_name in ds.data_vars:
    var = ds[var_name]
    print(f"\n🔹 Variable: {var_name}")
    print(f"   Dimensions: {var.dims}")
    print(f"   Shape: {var.shape}")
    print(f"   Attributes:")
    for attr_name, attr_value in var.attrs.items():
        print(f"     - {attr_name}: {attr_value}")

# Close the datasetds.close()
