import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import os
from geopy.distance import geodesic
from metpy.units import units
from metpy.calc import dewpoint_from_specific_humidity, equivalent_potential_temperature
from metpy.interpolate import cross_section

# -----------------------------
# 1. Input directories & transect
# -----------------------------
temp_dir = r'event\t_event.nc'   # ERA5 temperature NetCDF
q_dir = r'event\q_event.nc'      # ERA5 specific humidity NetCDF
output_dir = r'C:\Users\DJEV\Documents\Palawan_Shear_Line\ept_plot_output'
os.makedirs(output_dir, exist_ok=True)

# Transect start/end (lat, lon)
start = (8.0, 115.0)
end = (13.0, 127.0)

# -----------------------------
# 2. Load datasets
# -----------------------------
ds_t = xr.open_dataset(temp_dir).metpy.parse_cf()
ds_q = xr.open_dataset(q_dir).metpy.parse_cf()
ds = xr.merge([ds_t, ds_q])

# -----------------------------
# 3. Compute θₑ along transect
# -----------------------------
def compute_theta_e_cross(ds, start, end):
    """Compute equivalent potential temperature along a transect."""
    cross = cross_section(ds, start, end).set_coords(('latitude','longitude'))

    # Prepare units
    p = cross['pressure_level'] * 100 * units.pascal  # hPa → Pa
    T = cross['t'] * units.kelvin
    q = cross['q'] * units.dimensionless

    # Dewpoint & θe
    Td = dewpoint_from_specific_humidity(pressure=p, specific_humidity=q, temperature=T)
    theta_e = equivalent_potential_temperature(pressure=p, temperature=T, dewpoint=Td)

    # Clip extreme values to avoid huge numbers
    theta_e = theta_e.clip(min=250*units.kelvin, max=400*units.kelvin)

    # Assign θe to cross-section
    cross['theta_e'] = theta_e
    cross['theta_e'].attrs['units'] = 'K'

    print(f"θe range along transect: {theta_e.min().values:.1f} K to {theta_e.max().values:.1f} K")

    # Compute distance along transect
    lat = cross['latitude'].values
    lon = cross['longitude'].values
    dist = [0.0]
    for i in range(1, len(lat)):
        dist.append(dist[-1] + geodesic((lat[i-1], lon[i-1]), (lat[i], lon[i])).km)
    cross = cross.assign_coords(distance=('index', dist))  # Use 'index' as the x-axis

    return cross

# -----------------------------
# 4. Plot transect θₑ
# -----------------------------
def plot_theta_e_cross(cross, output_dir, transect_name='transect'):
    theta = cross['theta_e']

    # Remove extra dims
    dims_to_squeeze = [d for d in ['valid_time', 'number'] if d in theta.dims]
    theta_2d = theta.isel({d:0 for d in dims_to_squeeze}).values

    X, Y = np.meshgrid(cross['index'], cross['pressure_level'])  # no transpose
    Z = theta_2d  # ensure shape (pressure_level, distance)

    # X-axis labels (lat, lon)
    n_labels = 10
    total_points = len(cross['index'])
    label_idx = np.linspace(0, total_points-1, n_labels, dtype=int)
    xtick_labels = [f"{cross['latitude'][i].values:.2f}, {cross['longitude'][i].values:.2f}" 
                    for i in label_idx]

    plt.figure(figsize=(14,6))
    cf = plt.contourf(X, Y, Z, cmap='viridis')
    plt.gca().invert_yaxis()
    plt.colorbar(cf, label="θe [K]")
    plt.xlabel('Latitude, Longitude')
    plt.ylabel('Pressure [hPa]')
    plt.xticks(cross['index'][label_idx], xtick_labels, rotation=45, ha='right')
    plt.title(f'Equivalent Potential Temperature along {transect_name}')
    plt.tight_layout()

    outfile = os.path.join(output_dir, f'ept_{transect_name}.png')
    plt.savefig(outfile, dpi=150)
    plt.close()
    print(f"Plot saved: {outfile}")

# -----------------------------
# 5. Plot θₑ mean along latitude
# -----------------------------
def plot_theta_e_lat_mean(ds, output_dir):
    lat_min, lat_max = 8.0, 13.0
    lon_min, lon_max = 115.0, 127.0
    ds_sub = ds.sel(latitude=slice(lat_max, lat_min), longitude=slice(lon_min, lon_max))

    p = ds_sub['pressure_level'] * 100 * units.pascal
    T = ds_sub['t'] * units.kelvin
    q = ds_sub['q'] * units.dimensionless
    Td = dewpoint_from_specific_humidity(pressure=p, specific_humidity=q, temperature=T)
    theta_e = equivalent_potential_temperature(pressure=p, temperature=T, dewpoint=Td)
    theta_e = theta_e.clip(min=250*units.kelvin, max=400*units.kelvin)

    # Remove extra dims
    dims_to_squeeze = [d for d in ['valid_time','number'] if d in theta_e.dims]
    theta_2d = theta_e.isel({d:0 for d in dims_to_squeeze})

    # Average along longitude
    theta_lat = theta_2d.mean(dim='longitude')

    X, Y = np.meshgrid(theta_lat['latitude'], theta_lat['pressure_level'])
    plt.figure(figsize=(10,6))
    cf = plt.contourf(X, Y, theta_lat.values, cmap='viridis')
    plt.gca().invert_yaxis()
    plt.colorbar(cf, label="θe [K]")
    plt.xlabel('Latitude')
    plt.ylabel('Pressure [hPa]')
    plt.title('Equivalent Potential Temperature - Mean along Longitude')
    plt.tight_layout()

    outfile = os.path.join(output_dir, 'ept_mean_latitude.png')
    plt.savefig(outfile, dpi=150)
    plt.close()
    print(f"Plot saved: {outfile}")

# -----------------------------
# 6. Plot θₑ mean along longitude
# -----------------------------
def plot_theta_e_lon_mean(ds, output_dir):
    lat_min, lat_max = 8.0, 13.0
    lon_min, lon_max = 115.0, 127.0
    ds_sub = ds.sel(latitude=slice(lat_max, lat_min), longitude=slice(lon_min, lon_max))

    p = ds_sub['pressure_level'] * 100 * units.pascal
    T = ds_sub['t'] * units.kelvin
    q = ds_sub['q'] * units.dimensionless
    Td = dewpoint_from_specific_humidity(pressure=p, specific_humidity=q, temperature=T)
    theta_e = equivalent_potential_temperature(pressure=p, temperature=T, dewpoint=Td)
    theta_e = theta_e.clip(min=250*units.kelvin, max=400*units.kelvin)

    # Remove extra dims
    dims_to_squeeze = [d for d in ['valid_time','number'] if d in theta_e.dims]
    theta_2d = theta_e.isel({d:0 for d in dims_to_squeeze})

    # Average along latitude
    theta_lon = theta_2d.mean(dim='latitude')

    X, Y = np.meshgrid(theta_lon['longitude'], theta_lon['pressure_level'])
    plt.figure(figsize=(10,6))
    cf = plt.contourf(X, Y, theta_lon.values, cmap='viridis')
    plt.gca().invert_yaxis()
    plt.colorbar(cf, label="θe [K]")
    plt.xlabel('Longitude')
    plt.ylabel('Pressure [hPa]')
    plt.title('Equivalent Potential Temperature - Mean along Latitude')
    plt.tight_layout()

    outfile = os.path.join(output_dir, 'ept_mean_longitude.png')
    plt.savefig(outfile, dpi=150)
    plt.close()
    print(f"Plot saved: {outfile}")

# -----------------------------
# 7. Execute plotting
# -----------------------------
cross = compute_theta_e_cross(ds, start, end)
plot_theta_e_cross(cross, output_dir, transect_name='8N115E_to_13N127E')
plot_theta_e_lat_mean(ds, output_dir)
plot_theta_e_lon_mean(ds, output_dir)
