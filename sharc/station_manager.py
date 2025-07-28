# -*- coding: utf-8 -*-
"""
Created on Fri Feb  3 15:29:48 2017

@author: edgar
"""

import numpy as np

from sharc.support.enumerations import StationType
from sharc.station import Station
from sharc.antenna.antenna import Antenna
from sharc.mask.spectral_mask import SpectralMask


class StationManager(object):
    """
    This is the base class that manages an array of stations that will be
    used during a simulation. It acts like a container that vectorizes the
    station properties to speed up calculations.
    """

    def __init__(self, n):
        self.num_stations = n
        self.x = np.empty(n)  # x coordinate
        self.y = np.empty(n)  # y coordinate
        self.z = np.empty(n)  # z coordinate (includes height above ground)
        self.azimuth = np.empty(n)
        self.elevation = np.empty(n)
        self.height = np.empty(n)  # station height above ground
        self.idx_orbit = np.empty(n)
        self.indoor = np.zeros(n, dtype=bool)
        self.active = np.ones(n, dtype=bool)
        self.tx_power = np.empty(n)
        self.rx_power = np.empty(n)
        self.rx_interference = np.empty(n)  # Rx interferece in dBW
        self.ext_interference = np.empty(n)  # External interferece in dBW
        self.antenna = np.empty(n, dtype=Antenna)
        self.bandwidth = np.empty(n)  # Bandwidth in MHz
        self.noise_figure = np.empty(n)
        self.noise_temperature = np.empty(n)
        self.thermal_noise = np.empty(n)
        self.total_interference = np.empty(n)
        self.pfd_external = np.empty(n)  # External PFD in dBW/m²/MHz
        self.pfd_external_aggregated = np.empty(n)  # Aggregated External PFD in dBW/m²/MHz
        self.snr = np.empty(n)
        self.sinr = np.empty(n)
        self.sinr_ext = np.empty(n)
        self.inr = np.empty(n)  # INR in dBm/MHz
        self.pfd = np.empty(n)  # Powerflux density in dBm/m^2
        self.spectral_mask = np.empty(n, dtype=SpectralMask)
        self.center_freq = np.empty(n)
        self.station_type = StationType.NONE
        self.is_space_station = False
        self.intersite_dist = 0.0

    def get_station_list(self, id=None) -> list:
        if (id is None):
            id = range(self.num_stations)
        station_list = list()
        for i in id:
            station_list.append(self.get_station(i))
        return station_list

    def get_station(self, id) -> Station:
        station = Station()
        station.id = id
        station.x = self.x[id]
        station.y = self.y[id]
        station.z = self.z[id]
        station.azimuth = self.azimuth[id]
        station.elevation = self.elevation[id]
        station.height = self.height[id]
        station.indoor = self.indoor[id]
        station.active = self.active[id]
        station.tx_power = self.tx_power[id]
        station.rx_power = self.rx_power[id]
        station.rx_interference = self.rx_interference[id]
        station.ext_interference = self.ext_interference[id]
        station.antenna = self.antenna[id]
        station.bandwidth = self.bandwidth[id]
        station.noise_figure = self.noise_figure[id]
        station.noise_temperature = self.noise_temperature[id]
        station.thermal_noise = self.thermal_noise[id]
        station.total_interference = self.total_interference[id]
        station.snr = self.snr[id]
        station.sinr = self.sinr[id]
        station.sinr_ext = self.sinr_ext[id]
        station.inr = self.inr[id]
        station.station_type = self.station_type
        return station

    def get_distance_to(self, station) -> np.array:
        """Calculates the 2D distance between stations

        Parameters
        ----------
        station : StationManger
            Station to which calculate the distance

        Returns
        -------
        np.array
            Distance between stations
        """
        distance = np.empty([self.num_stations, station.num_stations])
        for i in range(self.num_stations):
            distance[i] = np.sqrt(
                np.power(self.x[i] - station.x, 2) +
                np.power(self.y[i] - station.y, 2),
            )
        return distance

    def get_3d_distance_to(self, station) -> np.array:
        """Calculates the 3D distance between stations

        Parameters
        ----------
        station : StationManager
            Station to which calculate the distance

        Returns
        -------
        np.array
            3D Distance between stations
        """
        distance = np.empty([self.num_stations, station.num_stations])
        for i in range(self.num_stations):
            distance[i] = np.sqrt(
                np.power(self.x[i] - station.x, 2) +
                np.power(self.y[i] - station.y, 2) +
                np.power(self.z[i] - station.z, 2)
            )
        return distance

    def get_dist_angles_wrap_around(self, station) -> np.array:
        """
        Calcualtes distances and angles using the wrap around technique
        Parameters:
            station (StationManager): station to which calculate the distances
                and angles
        Returns:
            distance_2D (np.array): 2D distance between stations
            distance_3D (np.array): 3D distance between stations
            phi (np.array): azimuth of pointing vector to other stations
            theta (np.array): elevation of pointing vector to other stations
        """
        # Initialize variables
        distance_3D = np.empty([self.num_stations, station.num_stations])
        distance_2D = np.inf * np.ones_like(distance_3D)
        cluster_num = np.zeros_like(distance_3D, dtype=int)

        # Cluster coordinates
        cluster_x = np.array([
            station.x,
            station.x + 3.5 * self.intersite_dist,
            station.x - 0.5 * self.intersite_dist,
            station.x - 4.0 * self.intersite_dist,
            station.x - 3.5 * self.intersite_dist,
            station.x + 0.5 * self.intersite_dist,
            station.x + 4.0 * self.intersite_dist,
        ])

        cluster_y = np.array([
            station.y,
            station.y + 1.5 *
            np.sqrt(3.0) * self.intersite_dist,
            station.y + 2.5 *
            np.sqrt(3.0) * self.intersite_dist,
            station.y + 1.0 *
            np.sqrt(3.0) * self.intersite_dist,
            station.y - 1.5 *
            np.sqrt(3.0) * self.intersite_dist,
            station.y - 2.5 *
            np.sqrt(3.0) * self.intersite_dist,
            station.y - 1.0 * np.sqrt(3.0) * self.intersite_dist,
        ])

        # Calculate 2D distance
        temp_distance = np.zeros_like(distance_2D)
        for k, (x, y) in enumerate(zip(cluster_x, cluster_y)):
            temp_distance = np.sqrt(
                np.power(x - self.x[:, np.newaxis], 2) +
                np.power(y - self.y[:, np.newaxis], 2),
            )
            is_shorter = temp_distance < distance_2D
            distance_2D[is_shorter] = temp_distance[is_shorter]
            cluster_num[is_shorter] = k

        # Calculate 3D distance
        distance_3D = np.sqrt(
            np.power(distance_2D, 2) +
            np.power(station.height - self.height[:, np.newaxis], 2),
        )

        # Calcualte pointing vector
        point_vec_x = cluster_x[cluster_num, np.arange(station.num_stations)] \
            - self.x[:, np.newaxis]
        point_vec_y = cluster_y[cluster_num, np.arange(station.num_stations)] \
            - self.y[:, np.newaxis]
        point_vec_z = station.height - self.height[:, np.newaxis]

        phi = np.array(
            np.rad2deg(
                np.arctan2(
                    point_vec_y, point_vec_x,
                ),
            ), ndmin=2,
        )
        theta = np.rad2deg(np.arccos(point_vec_z / distance_3D))

        return distance_2D, distance_3D, phi, theta

    def get_elevation(self, station) -> np.array:
        """
        Calculates the elevation angle between stations. Can be used for
        IMT stations.

        TODO: this implementation is essentialy the same as the one from
              get_elevation_angle (free-space elevation angle), despite the
              different matrix dimentions. So, the methods should be merged
              in order to reuse the source code
        """

        elevation = np.empty([self.num_stations, station.num_stations])

        for i in range(self.num_stations):
            distance = np.sqrt(
                np.power(self.x[i] - station.x, 2) +
                np.power(self.y[i] - station.y, 2),
            )
            rel_z = station.z - self.z[i]
            elevation[i] = np.degrees(np.arctan2(rel_z, distance))

        return elevation

    def get_pointing_vector_to(self, station) -> tuple:
        """calculate the pointing vector (angles) w.r.t. the other station

        Parameters
        ----------
        station : StationManager
            The other station to calculate the pointing vector

        Returns
        -------
        tuple
            phi, theta (phi is calculated with respect to x counter-clock-wise and
            theta is calculated with respect to z counter-clock-wise)
        """

        point_vec_x = station.x - self.x[:, np.newaxis]
        point_vec_y = station.y - self.y[:, np.newaxis]
        point_vec_z = station.z - self.z[:, np.newaxis]

        dist = self.get_3d_distance_to(station)

        phi = np.array(
            np.rad2deg(
                np.arctan2(
                    point_vec_y, point_vec_x,
                ),
            ), ndmin=2,
        )
        theta = np.rad2deg(np.arccos(point_vec_z / dist))

        return phi, theta

    def get_off_axis_angle(self, station) -> np.array:
        """
        Calculates the off-axis angle between this station and the input station
        """
        Az, b = self.get_pointing_vector_to(station)
        Az0 = self.azimuth

        a = 90 - self.elevation[:, np.newaxis]
        C = Az0[:, np.newaxis] - Az

        cos_phi = np.cos(np.radians(a)) * np.cos(np.radians(b)) \
                    + np.sin(np.radians(a)) * np.sin(np.radians(b)) * np.cos(np.radians(C))
        phi = np.arccos(
            # imprecision may accumulate enough for numbers to be slightly out of arccos range
            np.clip(cos_phi, -1., 1.)
        )
        phi_deg = np.degrees(phi)

        return phi_deg

    def is_imt_station(self) -> bool:
        """Whether this station is IMT or not

        Parameters
        ----------
        sta : StationManager
            The station that we're testing.

        Returns
        -------
        bool
            Whether this station is IMT or not
        """
        if self.station_type is StationType.IMT_BS or self.station_type is StationType.IMT_UE:
            return True
        else:
            return False


def copy_active_stations(stations: StationManager) -> StationManager:
    """Copy only the active stations from a StationManager object

    Parameters
    ----------
    stations : StationManager
        StationManager object to copy

    Returns
    -------
    StationManager
        A new StationManager object with only the active stations
    """
    act_sta = StationManager(np.sum(stations.active))
    for idx, active_idx in enumerate(np.where(stations.active)[0]):
        act_sta.x[idx] = stations.x[active_idx]
        act_sta.y[idx] = stations.y[active_idx]
        act_sta.z[idx] = stations.z[active_idx]
        act_sta.azimuth[idx] = stations.azimuth[active_idx]
        act_sta.elevation[idx] = stations.elevation[active_idx]
        act_sta.height[idx] = stations.height[active_idx]
        act_sta.indoor[idx] = stations.indoor[active_idx]
        act_sta.active[idx] = stations.active[active_idx]
        act_sta.tx_power[idx] = stations.tx_power[active_idx]
        act_sta.rx_power[idx] = stations.rx_power[active_idx]
        act_sta.rx_interference[idx] = stations.rx_interference[active_idx]
        act_sta.ext_interference[idx] = stations.ext_interference[active_idx]
        act_sta.antenna[idx] = stations.antenna[active_idx]
        act_sta.bandwidth[idx] = stations.bandwidth[active_idx]
        act_sta.noise_figure[idx] = stations.noise_figure[active_idx]
        act_sta.noise_temperature[idx] = stations.noise_temperature[active_idx]
        act_sta.thermal_noise[idx] = stations.thermal_noise[active_idx]
        act_sta.total_interference[idx] = stations.total_interference[active_idx]
        act_sta.snr[idx] = stations.snr[active_idx]
        act_sta.sinr[idx] = stations.sinr[active_idx]
        act_sta.sinr_ext[idx] = stations.sinr_ext[active_idx]
        act_sta.inr[idx] = stations.inr[active_idx]
        act_sta.pfd[idx] = stations.pfd[active_idx]
        act_sta.spectral_mask = stations.spectral_mask
        act_sta.center_freq[idx] = stations.center_freq[active_idx]
        act_sta.station_type = stations.station_type
        act_sta.is_space_station = stations.is_space_station
        act_sta.intersite_dist = stations.intersite_dist
    return act_sta
