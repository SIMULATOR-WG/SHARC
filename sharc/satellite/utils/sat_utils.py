import numpy as np
from sharc.parameters.constants import EARTH_RADIUS

EARTH_RADIUS_KM = EARTH_RADIUS / 1000


class WGS84Defs:
    """Constants for the WGS84 ellipsoid model."""
    SEMI_MAJOR_AXIS = 6378137.0  # Semi-major axis (in meters)
    SEMI_MINOR_AXIS = 6356752.3  # Semi-major axis (in meters)
    ECCENTRICITY = 8.1819190842622e-2  # WGS84 ellipsoid eccentricity
    FLATTENING = 0.0033528106647474805
    FIRST_ECCENTRICITY_SQRD = 6.69437999014e-3


def ecef2lla(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> tuple:
    """Coverts ECEF cartesian coordinates to lat long in WSG84 CRS.

    Parameters
    ----------
    x : np.ndarray
        x coordintate in meters
    y : np.ndarray
        y coordintate in meters
    z : np.ndarray
        x coordintate in meters

    Returns
    -------
    tuple (lat, long, alt)
        lat long and altitude in WSG84 format
    """
    # Longitude calculation
    lon = np.arctan2(y, x)

    # Iteratively solve for latitude and altitude
    p = np.sqrt(np.power(x, 2) + np.power(y, 2))
    lat = np.arctan2(z, p * (1 - WGS84Defs.ECCENTRICITY**2))  # Initial estimate for latitude
    for _ in range(5):  # Iteratively improve the estimate
        N = WGS84Defs.SEMI_MAJOR_AXIS / np.sqrt(1 - WGS84Defs.ECCENTRICITY**2 * np.sin(lat)**2)
        alt = p / np.cos(lat) - N
        lat = np.arctan2(z, p * (1 - WGS84Defs.ECCENTRICITY**2 * (N / (N + alt))))

    # Convert latitude and longitude from radians to degrees
    lat = np.degrees(lat)
    lon = np.degrees(lon)
    return lat, lon, alt


def lla2ecef(lat: np.ndarray, lng: np.ndarray, alt: np.ndarray) -> tuple:
    """Converts from geodetic WSG84 to ECEF coordinates

    Parameters
    ----------
    lat : np.ndarray
        latitude in degrees
    lng : np.ndarray
        longitute in degrees
    alt : np.ndarray
        altitude in meters

    Returns
    -------
    tuple
        x, y and z coordinates
    """
    lat = np.deg2rad(lat)
    lng = np.deg2rad(lng)
    n_phi = WGS84Defs.SEMI_MAJOR_AXIS / np.sqrt(1 - WGS84Defs.FIRST_ECCENTRICITY_SQRD * np.sin(lat)**2)
    x = (n_phi + alt) * np.cos(lat) * np.cos(lng)
    y = (n_phi + alt) * np.cos(lat) * np.sin(lng)
    z = ((1 - WGS84Defs.FLATTENING)**2 * n_phi + alt) * np.sin(lat)

    return x, y, z


def calc_elevation(Le: np.ndarray,
                   Ls: np.ndarray,
                   le: np.ndarray,
                   ls: np.ndarray,
                   sat_height: np.ndarray) -> np.ndarray:
    """Calculates the elevation angle from the earth station
    to space station, given earth and space station coordinates.
    Negative elevation angles means the space stations is not visible from Earth station.

    Parameters
    ----------
    Le : (ndarray)
        latitudes of the earth station
    Ls : (ndarray)
        latitudes of the space station
    le : (ndarray)
        longitudes of the earth station
    ls : (ndarray)
        latitudes of the space station
    sat_height : (ndarray)
        space station altitudes

    Returns
    -------
    (ndarray)
        array of elevation angles from the earth station in degrees.
    """
    Le = np.radians(Le)
    Ls = np.radians(Ls)
    le = np.radians(le)
    ls = np.radians(ls)
    gamma = np.arccos(
        np.cos(Le) * np.cos(Ls) * np.cos(ls - le) + np.sin(Le) * np.sin(Ls)
    )
    rs = EARTH_RADIUS_KM + sat_height
    slant = np.sqrt(rs**2 + EARTH_RADIUS_KM**2 - 2 * rs * EARTH_RADIUS_KM * np.cos(gamma))
    elev_angle = np.arccos((slant**2 + EARTH_RADIUS_KM**2 - rs**2) / \
                           (2 * slant * EARTH_RADIUS_KM)) - np.pi / 2

    return np.degrees(elev_angle)


if __name__ == "__main__":
    r1 = ecef2lla(7792.1450, 0, 0)
    print(r1)
    r2 = lla2ecef(r1[0], r1[1], r1[2])
    print(r2)
