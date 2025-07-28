import numpy as np
import pandas as pd
from scipy.special import jv
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import time

# START OF TIMER
start_time = time.time()

# CONSTANTS
EARTH_RADIUS_KM = 6378.145  # radius of the Earth, in km
KEPLER_CONST = 398601.8  # Kepler's constant, in km^3/s^2
EARTH_ROTATION_RATE = 2 * np.pi / (24 * 3600)  # earth's average rotation rate, in rad/s

# INPUTS - NGSO SYSTEM

# Globalstar
Nsp = 6  # number of satellites in the orbital plane (A.4.b.4.b)
Np = 8  # number of orbital planes (A.4.b.2)
phasing = 7.5  # satellite phasing between planes, in degrees
long_asc = 0  # initial longitude of ascending node of the first plane, in degrees (A.4.b.4.j)
omega = 0  # argument of perigee, in degrees (A.4.b.4.i)
delta = 52  # orbital plane inclination, in degrees (A.4.b.4.a)
hp = 1414  # altitude of perigee in km (A.4.b.4.e)
ha = 1414  # altitude of apogee in km (A.4.b.4.d)
Mo = 0  # initial mean anomaly for first satellite of first plane, in degrees

# # Molniya
# Nsp = 1  # number of satellites in the orbital plane (A.4.b.4.b)
# Np = 9  # number of orbital planes (A.4.b.2)
# phasing = 80 # Satellite phasing between planes (degrees)
# long_asc = 90  # longitude of ascending node in degrees (A.4.b.4.j) or RAAN
# omega = -90  # argument of perigee in degrees
# delta = 63.4  # inclination of the orbital plane in degrees
# hp = 950  # altitude at perigee in km
# ha = 39520  # altitude at apogee in km
# Mo = 0  # initial mean anomaly for first satellite of first plane

a = (hp + ha + 2 * EARTH_RADIUS_KM) / 2  # semi-major axis, in km
e = (ha + EARTH_RADIUS_KM - a) / a  # orbital eccentricity (e)
P = 2 * np.pi * np.sqrt(a ** 3 / KEPLER_CONST)  # orbital period, in seconds
beta = 360 / Nsp  # satellite separation angle in the plane (degrees)
psi = 360 / Np  # angle between plane intersections with the equatorial plane (degrees)

# Initial mean anomalies for all the satellites
M_o = (Mo + np.arange(Nsp) * beta + np.arange(Np)[:, None] * phasing) % 360  # shape (Np, Nsp)
M0 = np.radians(M_o.flatten())  # shape (Np*Nsp,)

# Initial longitudes of ascending node for all the planes
Omega_o = (long_asc + np.arange(Nsp) * 0 + np.arange(Np)[:, None] * psi) % 360  # shape (Np, Nsp)
Omega0 = np.radians(Omega_o.flatten())  # shape (Np*Nsp,)

# INPUTS - EARTH SYSTEM

# station position vector, in km
lat_es = [0, -10]  # latitudes of stations, in degrees
lon_es = [0, -30]  # longitudes of stations, in degrees
alt_es = [0, 0]  # altitudes of stations, in meters
rx = EARTH_RADIUS_KM + alt_es[0] * 10**-3
px = rx * np.cos(np.radians(lat_es[0])) * np.cos(np.radians(lon_es[0]))
py = rx * np.cos(np.radians(lat_es[0])) * np.sin(np.radians(lon_es[0]))
pz = rx * np.sin(np.radians(lat_es[0]))
p = np.array([px, py, pz])

# INPUT - TIME OF SIMULATION

interval = 5  # time interval, in seconds
nP = 4  # number of orbital periods of simulation
t = np.arange(0, nP * P + interval, interval)  # vector of time, in seconds

# CUSTOM FUNCTIONS


def wrap2pi(angle_rad):
    """
    Adjusts an angle in radians to the interval (-π, π).
    """
    return (angle_rad + np.pi) % (2 * np.pi) - np.pi


def eccentric_anomaly(e, M, terms=40, mod_2pi=True):
    """
    Calculate the eccentric anomaly E given the eccentricity e and
    mean anomaly M using a Bessel series expansion.

    Parameters:
    ----------
    e: float
        Eccentricity of the orbit (0 <= e < 1).
    M: np.ndarray
        Mean anomaly in radians; can be any shape.
    terms: int, optional
        Number of terms for the Bessel series expansion (default is 40).
    mod_2pi: bool, optional
        Whether to return E modulo 2π (default is True).

    Returns:
    -------
    np.ndarray
        Eccentric anomaly E in radians, with the same shape as M.
    """
    # Handle edge case where e is near zero, return M directly.
    if np.isclose(e, 0):
        return M

    # Prepare for Bessel series expansion
    n = np.arange(1, terms + 1)[:, None, None]  # Shape (terms, 1, 1)
    M_expanded = M[None, ...]  # Add a new dimension to M, for broadcasting

    # Calculate series sum using Bessel functions and sine terms
    series_sum = np.sum((1 / n) * jv(n, n * e) * np.sin(n * M_expanded), axis=0)

    # Calculate eccentric anomaly
    E = M + 2 * series_sum

    # Apply modulo operation if specified
    if mod_2pi:
        E = np.mod(E, 2 * np.pi)

    return E


def keplerian2eci(a, e, delta, Omega, omega, nu):
    """
    Calculate the position vector in ECI coordinates for each RAAN and true anomaly.

    Parameters:
    ----------
    a : float
        Semi-major axis in km.
    e : float
        Eccentricity (0 <= e < 1).
    delta : float
        Orbital inclination in degrees.
    Omega : np.ndarray
        Right ascension of the ascending node (RAAN) in degrees. Shape = (N,)
    omega : float
        Argument of perigee in degrees.
    nu : np.ndarray
        True anomaly in degrees. Shape = (N, length(t))

    Returns:
    -------
    r_eci : np.ndarray
        Position vector in ECI coordinates. Shape = (3, N, length(t))
    """
    # Convert angles from degrees to radians
    delta_rad = np.radians(delta)
    Omega_rad = np.radians(Omega)
    omega_rad = np.radians(omega)
    nu_rad = np.radians(nu)

    # Compute gamma (angle between satellite position and ascending node in the orbital plane)
    gamma = (nu_rad + omega_rad) % (2 * np.pi)

    # Trigonometric calculations
    cos_gamma = np.cos(gamma)
    sin_gamma = np.sin(gamma)
    cos_raan = np.cos(-Omega_rad)[:, np.newaxis]  # Shape (N, 1) for broadcasting
    sin_raan = np.sin(-Omega_rad)[:, np.newaxis]  # Shape (N, 1) for broadcasting
    cos_incl = np.cos(delta_rad)
    sin_incl = np.sin(delta_rad)

    # Calculate radius for each true anomaly
    r = a * (1 - e ** 2) / (1 + e * np.cos(nu_rad))

    # Position in ECI coordinates
    x = r * (cos_gamma * cos_raan - sin_gamma * sin_raan * cos_incl)
    y = r * (cos_gamma * sin_raan + sin_gamma * cos_raan * cos_incl)
    z = r * sin_gamma * sin_incl

    # Stack to form the ECI position vector with shape (3, N, length(t))
    r_eci = np.array([x, y, z])

    return r_eci


def eci2ecef(t, r_eci):
    """
    Convert coordinates from Earth-Centered Inertial (ECI) to Earth-Centered Earth-Fixed (ECEF)
    for each time step, accounting for Earth's rotation.

    Parameters:
    ----------
    t : np.ndarray
        A 1D array of time instants in seconds since the reference epoch. Shape = (T,)
    r_eci : np.ndarray
        A 3D array representing the position vectors in ECI coordinates (km). Shape = (3, N, T)

    Returns:
    -------
    r_ecef : np.ndarray
        A 3D array representing the position vectors in ECEF coordinates (km). Shape = (3, N, T)
    """
    # Earth's rotation rate in radians per second (WGS-84 value)
    ωe = 2 * np.pi / (24 * 3600)  # rad/s

    # Calculate the rotation angle θ for each time instant in `t`
    θ = wrap2pi(ωe * t)  # Shape (T,)

    # Create the rotation matrices for each θ in the array
    cos_θ = np.cos(θ)
    sin_θ = np.sin(θ)

    # Rotation matrices for each time step `t`, shape (T, 3, 3)
    S = np.array([
        [cos_θ, sin_θ, np.zeros_like(θ)],
        [-sin_θ, cos_θ, np.zeros_like(θ)],
        [np.zeros_like(θ), np.zeros_like(θ), np.ones_like(θ)]
    ]).transpose(2, 0, 1)  # Shape (T, 3, 3)

    # Transpose `r_eci` for easier matrix multiplication, shape (T, 3, N)
    r_eci_t = r_eci.transpose(2, 0, 1)

    # Apply the rotation for each time step using `einsum`
    r_ecef_t = np.einsum('tij,tjk->tik', S, r_eci_t)  # Shape (T, 3, N)

    # Transpose back to shape (3, N, T) for the final result
    r_ecef = r_ecef_t.transpose(1, 2, 0)

    return r_ecef


def plot_ground_tracks(theta_deg, phi_deg, planes=None, satellites=None, title="Satellite Ground Tracks"):
    """
    Plots the satellite ground tracks as points based on latitude and longitude data.

    Parameters:
        theta_deg (ndarray): Array of satellite latitudes in degrees, shape (num_planes * num_satellites_per_plane,
        num_timesteps).
        phi_deg (ndarray): Array of satellite longitudes in degrees, shape (num_planes * num_satellites_per_plane,
        num_timesteps).
        planes (list, optional): List of plane indices to plot. If None, all planes are plotted.
        satellites (list, optional): List of satellite indices to plot. If None, all satellites are plotted.
        title (str): Title for the plot.
    """
    # Determine the number of planes and satellites in the data
    num_planes = len(np.unique(planes)) if planes is not None else theta_deg.shape[0]
    num_satellites_per_plane = theta_deg.shape[0] // num_planes

    # Set up the plot with Cartopy
    fig, ax = plt.subplots(figsize=(12, 8), subplot_kw={'projection': ccrs.PlateCarree()})
    ax.set_global()
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    ax.gridlines(draw_labels=True)

    # Loop through each plane and satellite, plotting ground tracks as points
    for plane_idx in range(num_planes):
        # Skip plane if it's not selected
        if planes and plane_idx + 1 not in planes:
            continue

        for sat_idx in range(num_satellites_per_plane):
            # Calculate satellite index based on plane and satellite number
            satellite_global_idx = plane_idx * num_satellites_per_plane + sat_idx

            # Skip satellite if it's not selected
            if satellites and sat_idx + 1 not in satellites:
                continue

            # Plot the ground track points for the selected satellite
            ax.scatter(
                phi_deg[satellite_global_idx, :], theta_deg[satellite_global_idx, :],
                label=f'Plane {plane_idx + 1}, Satellite {sat_idx + 1}', s=1,  # Size of each point
                transform=ccrs.PlateCarree()
            )

    # Add legend and title
    ax.set_title(title)
    ax.legend(loc='upper right', markerscale=5, fontsize='small')
    plt.show()


if __name__ == "__main__":
    # POSITION CALCULATION IN THE ORBITAL PLANE (v and r)

    # Mean anomaly (M)
    mean_anomaly = (M0[:, None] + (2 * np.pi / P) * t) % (2 * np.pi)

    # Eccentric anomaly (E)
    eccentric_anom = eccentric_anomaly(e, mean_anomaly)

    # True anomaly (v)
    v = 2 * np.arctan(np.sqrt((1 + e) / (1 - e)) * np.tan(eccentric_anom / 2))
    v = np.mod(v, 2 * np.pi)

    # Distance of the satellite to Earth's center (r)
    r = a * (1 - e ** 2) / (1 + e * np.cos(v))

    # True anomaly relative to the line of nodes (gamma)
    gamma = wrap2pi(v + np.radians(omega))  # gamma in the interval [-pi, pi]

    # Latitudes of the satellites, in radians (theta)
    theta = np.arcsin(np.sin(gamma) * np.sin(np.radians(delta)))

    # Longitude variation due to angular displacement, in radians (phiS)
    phiS = np.arccos(np.cos(gamma) / np.cos(theta)) * np.sign(gamma)

    # Longitudes of the ascending node (OmegaG)
    OmegaG = (Omega0[:, None] + EARTH_ROTATION_RATE * t)  # shape (Np*Nsp, len(t))
    OmegaG = wrap2pi(OmegaG)

    # POSITION CALCULATION IN ECEF COORDINATES - ITU-R S.1503

    r_eci = keplerian2eci(a, e, delta, np.degrees(Omega0), omega, np.degrees(v))
    r_ecef = eci2ecef(t, r_eci)
    sx, sy, sz = r_ecef[0], r_ecef[1], r_ecef[2]
    lat = np.degrees(np.arcsin(sz / r))
    lon = np.degrees(np.arctan2(sy, sx))

    # MAP
    # static map
    plot_ground_tracks(lat, lon, planes=[1, 2, 3, 4, 5, 6, 7, 8], satellites=[1])

    # DATAFRAME WITH INITIAL DATA
    data = []
    for i, time_step in enumerate(t):
        for plane in range(Np):
            for satellite in range(Nsp):
                sat_index = plane * Nsp + satellite  # Calculate the index for flattened arrays
                # Append data for each satellite at each time step
                data.append([
                    time_step,  # Time
                    plane + 1,  # Plane number
                    satellite + 1,  # Satellite number
                    M_o[plane, satellite],  # Mean anomaly (M_o) in degrees
                    Omega_o[plane, satellite],  # Longitude of ascending node (Omega_o) in degrees
                    r[sat_index, i],  # distance in km
                    lat[sat_index, i],  # Latitude in degrees
                    lon[sat_index, i],  # Longitude in degrees
                    sx[sat_index, i],  # Cartesian x position
                    sy[sat_index, i],  # Cartesian y position
                    sz[sat_index, i]  # Cartesian z position
                ])

    df = pd.DataFrame(data, columns=['time',
                                     'plane',
                                     'satellite',
                                     'M_o',
                                     'Omega_o',
                                     'r',
                                     'lat',
                                     'lon',
                                     'sx',
                                     'sy',
                                     'sz'])

    # END OF TIMER
    end_time = time.time()

    # Calculate elapsed time
    elapsed_time = end_time - start_time
    print(f"Execution time: {elapsed_time:.4f} seconds")
