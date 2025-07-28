import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from scipy.special import jv


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
    e : float
        Eccentricity of the orbit (0 <= e < 1).
    M : np.ndarray
        Mean anomaly in radians; can be any shape.
    terms : int, optional
        Number of terms for the Bessel series expansion (default is 40).
    mod_2pi : bool, optional
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
