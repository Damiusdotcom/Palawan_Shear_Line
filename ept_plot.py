import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import os
from geopy.distance import geodesic

from metpy.interpolate import cross_section
from metpy.units import units
from metpy.calc import dewpoint_from_specific_humidity, equivalent_potential_temperature

# -----------------------------
# 1. Input directories & transect
# -----------------------------
temp_dir = r'event\t_event.nc'   # ERA5 temperature NetCDF
q_dir = r'event\q_event.nc'               # ERA5 specific humidity NetCDF

output_dir = r'C:\Users\DJEV\Documents\Palawan_Shear_Line\ept_plot_output'
os.makedirs(output_dir, exist_ok=True)

# Transect start/end (lat, lon)
start = (8.0, 115.0)   # 8N, 115E
end = (13.0, 127.0)    # 13N, 127E

# -----------------------------
# 2. Load datasets
# -----------------------------
ds_t = xr.open_dataset(temp_dir).metpy.parse_cf()
ds_q = xr.open_dataset(q_dir).metpy.parse_cf()

# Merge datasets into one
ds = xr.merge([ds_t, ds_q])

# -----------------------------
# 3. Compute θₑ along transect
# -----------------------------
def compute_theta_e_cross(ds, start, end):
    """
    Computes equivalent potential temperature along a transect.
    Returns cross-section DataArray with theta_e.
    """
    # Extract cross-section
    cross = cross_section(ds, start, end).set_coords(('latitude', 'longitude'))

    # Variables with proper units
    p_full = xr.ones_like(cross['t']) * (cross['pressure_level'] * 100 * units.pascal)  # hPa → Pa
    T = cross['t'] * units.kelvin
    q = cross['q'] * units.dimensionless  # kg/kg

    # Compute dewpoint from specific humidity
    Td = dewpoint_from_specific_humidity(pressure=p_full, specific_humidity=q, temperature=T)

    # Compute equivalent potential temperature
    theta_e = equivalent_potential_temperature(p_full, T, Td)

    # Attach to cross-section
    cross['theta_e'] = theta_e
    cross['theta_e'].attrs['units'] = 'K'

    # Compute distance along transect (km)
    lat = cross['latitude'].values
    lon = cross['longitude'].values
    dist = [0.0]
    for i in range(1, len(lat)):
        dist.append(dist[-1] + geodesic((lat[i-1], lon[i-1]), (lat[i], lon[i])).km)
    cross = cross.assign_coords(distance=('points', dist))

    return cross

cross = compute_theta_e_cross(ds, start, end)

# -----------------------------
# 4. Plotting function
# -----------------------------
def plot_theta_e_cross(cross, output_dir, transect_name='transect'):
    # Select 2D theta_e
    theta_e_2d = cross['theta_e'].squeeze()
    if 'valid_time' in theta_e_2d.dims:
        theta_e_2d = theta_e_2d.isel(valid_time=0)
    if 'number' in theta_e_2d.dims:
        theta_e_2d = theta_e_2d.isel(number=0)

    # Use points as X-axis for contourf
    points = np.arange(len(cross['latitude']))
    X, Y = np.meshgrid(points, cross['pressure_level'])

    # Create X-axis labels for lat/lon every nth point
    n_labels = 10
    xticks = np.linspace(0, len(points)-1, n_labels, dtype=int)
    xtick_labels = [f"({cross['latitude'][i].values:.2f}, {cross['longitude'][i].values:.2f})" for i in xticks]

    # Create Y-axis labels using actual pressure levels
    ytick_labels = [f"{int(p)}" for p in cross['pressure_level'].values]

    # Plot
    plt.figure(figsize=(14,6))
    plt.contourf(X, Y, theta_e_2d, cmap='viridis')
    plt.gca().invert_yaxis()  # Pressure decreases upwards
    plt.colorbar(label=f"Equivalent Potential Temperature [{theta_e_2d.units}]")
    plt.xlabel('Latitude, Longitude')
    plt.ylabel('Pressure [hPa]')
    plt.xticks(xticks, xtick_labels, rotation=45, ha='right')
    plt.yticks(cross['pressure_level'].values, ytick_labels)
    plt.title(f'Equivalent Potential Temperature along {transect_name}')
    plt.tight_layout()

    # Save PNG
    outfile = os.path.join(output_dir, f'ept_{transect_name}.png')
    plt.savefig(outfile, dpi=150)
    plt.close()
    print(f"Plot saved: {outfile}")


plot_theta_e_cross(cross, output_dir, transect_name='8N115E_to_13N127E')
