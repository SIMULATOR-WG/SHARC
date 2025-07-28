"""Implements a Space Station Orbit model as described in Rec. ITU-R S.1325-3
"""

import numpy as np

from sharc.satellite.ngso.custom_functions import wrap2pi, eccentric_anomaly, keplerian2eci, eci2ecef, plot_ground_tracks
from sharc.satellite.utils.sat_utils import ecef2lla
from sharc.satellite.ngso.constants import EARTH_RADIUS_KM, KEPLER_CONST, EARTH_ROTATION_RATE


class OrbitModel():
    """Orbit Model for satellite positions."""

    def __init__(self,
                 Nsp: int,
                 Np: int,
                 phasing: float,
                 long_asc: float,
                 omega: float,
                 delta: float,
                 hp: float,
                 ha: float,
                 Mo: float):
        """Instantiates and OrbitModel object from the Orbit parameters as specified in S.1529.

        Parameters
        ----------
        Nsp : int
            number of satellites in the orbital plane (A.4.b.4.b)
        Np : int
            number of orbital planes (A.4.b.2)
        phasing : float
            satellite phasing between planes, in degrees
        long_asc : float
            initial longitude of ascending node of the first plane, in degrees
        omega : float
            argument of perigee, in degrees
        delta : float
            orbital plane inclination, in degrees
        perigee_alt_km : float
            altitude of perigee in km
        ha : float
            altitude of apogee in km
        Mo : float
            initial mean anomaly for first satellite of first plane, in degrees
        """
        self.Nsp = Nsp
        self.Np = Np
        self.phasing = phasing
        self.long_asc = long_asc
        self.omega = omega
        self.delta = delta
        self.perigee_alt_km = hp
        self.apogee_alt_km = ha
        self.Mo = Mo

        # Derive other orbit parameters
        self.semi_major_axis = (hp + ha + 2 * EARTH_RADIUS_KM) / 2  # semi-major axis, in km
        self.eccentricity = (ha + EARTH_RADIUS_KM - self.semi_major_axis) / self.semi_major_axis  # orbital eccentricity (e)
        self.orbital_period_sec = 2 * np.pi * np.sqrt(self.semi_major_axis ** 3 / KEPLER_CONST)  # orbital period, in seconds
        self.sat_sep_angle_deg = 360 / self.Nsp  # satellite separation angle in the plane (degrees)
        self.orbital_plane_inclination = 360 / self.Np  # angle between plane intersections with the equatorial plane (degrees)

        # Initial mean anomalies for all the satellites
        initial_mean_anomalies_deg = (Mo + np.arange(Nsp) * self.sat_sep_angle_deg + np.arange(self.Np)[:, np.newaxis] *
                                      self.phasing) % 360
        self.initial_mean_anomalies_rad = np.radians(initial_mean_anomalies_deg.flatten())

        # Initial longitudes of ascending node for all the planes
        inital_raan = (self.long_asc * np.ones(self.Nsp) + np.arange(self.Np)[:, np.newaxis] *
                       self.orbital_plane_inclination) % 360
        self.inital_raan_rad = np.radians(inital_raan.flatten())

        # Calculated variables
        self.mean_anomaly = None  # computed mean anomalies
        self.raan_rad = None  # computed longitudes of the ascending node
        self.eccentric_anomaly = None  # computed eccentric anomalies
        self.true_anomaly = None  # computed true anomalies
        self.distance = None  # computed distance to Earth's center
        self.true_anomaly_rel_line_nodes = None  # computed true anomaly relative to the line of nodes

    def get_satellite_positions_time_interval(self, initial_time_secs=0, interval_secs=5, n_periods=4) -> dict:
        """
        Return the orbit positions vector.

        Parameters
        ----------
        initial_time_secs : int, optional
            initial time instant in seconds, by default 0
        interval_secs : int, optional
            time interval between points, by default 5
        n_periods : int, optional
            number of orbital peridos, by default 4

        Returns
        -------
        dict
            A dictionary with satellite positions in spherical and ecef coordinates.
                lat, lon, sx, sy, sz
        """
        t = np.arange(initial_time_secs, n_periods * self.orbital_period_sec + interval_secs, interval_secs)
        return self.__get_satellite_positions(t)

    def get_orbit_positions_time_instant(self, time_instant_secs=0) -> dict:
        """Returns satellite positions in a determined time instant in seconds"""
        t = np.array([time_instant_secs])
        return self.__get_satellite_positions(t)

    def get_orbit_positions_random(self, rng: np.random.RandomState, n_samples=1) -> dict:
        """Returns satellite positions in a random time instant in seconds.
                Parameters
                ----------
                rng : np.random.RandomState
                    Random number generator for reproducibility
                n_samples : int
                    Number of random samples to generate, by default 1
                Returns
                -------
                dict
                    A dictionary with satellite positions in spherical and ecef coordinates.
                        lat, lon, sx, sy, sz
        """

        # return self.__get_satellite_positions(rng.random_sample(1) * 1000 * self.orbital_period_sec)
        # Mean anomaly (M)
        self.mean_anomaly = (self.initial_mean_anomalies_rad[:, None] +
                             2 * np.pi * rng.random_sample(n_samples)) % (2 * np.pi)

        # Eccentric anomaly (E)
        self.eccentric_anom = eccentric_anomaly(self.eccentricity, self.mean_anomaly)

        # True anomaly (v)
        self.true_anomaly = 2 * np.arctan(np.sqrt((1 + self.eccentricity) / (1 - self.eccentricity)) *
                                          np.tan(self.eccentric_anom / 2))

        self.true_anomaly = np.mod(self.true_anomaly, 2 * np.pi)

        # Distance of the satellite to Earth's center (r)
        r = self.semi_major_axis * (1 - self.eccentricity ** 2) / (1 + self.eccentricity * np.cos(self.true_anomaly))

        # True anomaly relative to the line of nodes (gamma)
        self.true_anomaly_rel_line_nodes = wrap2pi(self.true_anomaly + np.radians(self.omega))  # gamma in the interval [-pi, pi]

        # Latitudes of the satellites, in radians (theta)
        # theta = np.arcsin(np.sin(gamma) * np.sin(np.radians(self.delta)))

        # Longitude variation due to angular displacement, in radians (phiS)
        # phiS = np.arccos(np.cos(gamma) / np.cos(theta)) * np.sign(gamma)

        # Longitudes of the ascending node (OmegaG)
        raan_rad = (self.inital_raan_rad[:, None] + 2 * np.pi * rng.random_sample(n_samples))  # shape (Np*Nsp, len(t))
        raan_rad = wrap2pi(raan_rad)

        # POSITION CALCULATION IN ECEF COORDINATES - ITU-R S.1503
        r_eci = keplerian2eci(self.semi_major_axis,
                              self.eccentricity,
                              self.delta,
                              np.degrees(self.inital_raan_rad),
                              self.omega,
                              np.degrees(self.true_anomaly))

        r_ecef = eci2ecef(2 * np.pi * rng.random_sample(n_samples) / EARTH_ROTATION_RATE, r_eci)
        sx, sy, sz = r_ecef[0], r_ecef[1], r_ecef[2]
        lat = np.degrees(np.arcsin(sz / r))
        lon = np.degrees(np.arctan2(sy, sx))
        # (lat, lon, _) = ecef2lla(sx, sy, sz)

        pos_vector = {
            'lat': lat,
            'lon': lon,
            'sx': sx,
            'sy': sy,
            'sz': sz
        }
        return pos_vector

    def __get_satellite_positions(self, t: np.array) -> dict:
        """Returns the Satellite positins (both lat long and ecef) for a given time vector within the orbit period.

        Parameters
        ----------
        t : np.array
            time instants inside the orbit period in seconds

        Returns
        -------
        dict
            A dictionary with satellite positions in spherical and ecef coordinates.
                lat, lon, sx, sy, sz
        """
        # Mean anomaly (M)
        self.mean_anomaly = (self.initial_mean_anomalies_rad[:, None] +
                             (2 * np.pi / self.orbital_period_sec) * t) % (2 * np.pi)

        # Eccentric anomaly (E)
        self.eccentric_anom = eccentric_anomaly(self.eccentricity, self.mean_anomaly)

        # True anomaly (v)
        self.true_anomaly = 2 * np.arctan(np.sqrt((1 + self.eccentricity) / (1 - self.eccentricity)) *
                                          np.tan(self.eccentric_anom / 2))
        self.true_anomaly = np.mod(self.true_anomaly, 2 * np.pi)

        # Distance of the satellite to Earth's center (r)
        r = self.semi_major_axis * (1 - self.eccentricity ** 2) / (1 + self.eccentricity * np.cos(self.true_anomaly))

        # True anomaly relative to the line of nodes (gamma)
        self.true_anomaly_rel_line_nodes = wrap2pi(self.true_anomaly + np.radians(self.omega))  # gamma in the interval [-pi, pi]

        # Latitudes of the satellites, in radians (theta)
        # theta = np.arcsin(np.sin(gamma) * np.sin(np.radians(self.delta)))

        # Longitude variation due to angular displacement, in radians (phiS)
        # phiS = np.arccos(np.cos(gamma) / np.cos(theta)) * np.sign(gamma)

        # Longitudes of the ascending node (OmegaG)
        raan_rad = (self.inital_raan_rad[:, None] + EARTH_ROTATION_RATE * t)  # shape (Np*Nsp, len(t))
        raan_rad = wrap2pi(raan_rad)

        # POSITION CALCULATION IN ECEF COORDINATES - ITU-R S.1503
        r_eci = keplerian2eci(self.semi_major_axis,
                              self.eccentricity,
                              self.delta,
                              np.degrees(self.inital_raan_rad),
                              self.omega,
                              np.degrees(self.true_anomaly))

        r_ecef = eci2ecef(t, r_eci)
        sx, sy, sz = r_ecef[0], r_ecef[1], r_ecef[2]
        lat = np.degrees(np.arcsin(sz / r))
        lon = np.degrees(np.arctan2(sy, sx))
        # (lat, lon, _) = ecef2lla(sx, sy, sz)

        pos_vector = {
            'lat': lat,
            'lon': lon,
            'sx': sx,
            'sy': sy,
            'sz': sz
        }
        return pos_vector


def main():
    """Main function to test the OrbitModel class and plot ground tracks.

    This function creates an instance of the OrbitModel class with specified parameters,
    retrieves satellite positions over a specified time interval, and plots the ground tracks
    of the satellites using Plotly.
    """
    import plotly.graph_objects as go

    orbit_params = {
        "Nsp": 1,
        "Np": 28,
        "phasing": 1.5,
        "long_asc": 0,
        "omega": 0,
        "delta": 53,
        "hp": 525,
        "ha": 525,
        "Mo": 0
    }

    print("Orbit parameters:")
    print(orbit_params)

    # Instantiate the OrbitModel
    orbit_model = OrbitModel(**orbit_params)

    # Get satellite positions over time
    positions = orbit_model.get_satellite_positions_time_interval(n_periods=1)

    # Extract latitude and longitude
    latitudes = positions['lat']
    longitudes = positions['lon']

    # Create a plotly figure
    fig = go.Figure()

    # Add traces for each satellite
    for i in range(latitudes.shape[0]):
    # for i in range(10, 20):
        fig.add_trace(go.Scattergeo(
            lon=longitudes[i],
            lat=latitudes[i],
            mode='lines',
            name=f'Satellite {i + 1}'
        ))

    # Update layout for better visualization
    fig.update_layout(
        title="Satellite Ground Tracks",
        showlegend=False,
        geo=dict(
            projection_type="equirectangular",
            showland=True,
            landcolor="rgb(243, 243, 243)",
            countrycolor="rgb(204, 204, 204)"
        )
    )

    # Show the plot
    fig.show()


if __name__ == "__main__":
    main()
