import xarray as xr
import os

# -----------------------------
# 1. Input files
# -----------------------------
temp_file = r'event\t_event.nc'
q_file = r'event\q_event.nc'

# -----------------------------
# 2. Load datasets
# -----------------------------
ds_t = xr.open_dataset(temp_file)
ds_q = xr.open_dataset(q_file)

# -----------------------------
# 3. Diagnostic function
# -----------------------------
def diagnose_dataset(ds, name='Dataset'):
    print(f"\n==== {name} ====")
    print("Variables:")
    for var in ds.data_vars:
        print(f"  {var} - dims: {ds[var].dims}, shape: {ds[var].shape}")
        print(f"    Attributes: {ds[var].attrs}")
    print("\nCoordinates:")
    for coord in ds.coords:
        print(f"  {coord} - dims: {ds[coord].dims}, shape: {ds[coord].shape}")
        print(f"    Attributes: {ds[coord].attrs}")
    print("\nGlobal attributes:")
    for attr in ds.attrs:
        print(f"  {attr}: {ds.attrs[attr]}")

# -----------------------------
# 4. Run diagnostics
# -----------------------------
diagnose_dataset(ds_t, name='Temperature Dataset (t)')
diagnose_dataset(ds_q, name='Specific Humidity Dataset (q)')
