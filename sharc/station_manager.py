# -*- coding: utf-8 -*-
"""
Created on Fri Feb  3 15:29:48 2017

@author: edgar
"""

import numpy as np
import math

from sharc.support.enumerations import StationType
from sharc.station import Station
from sharc.antenna.antenna import Antenna

class StationManager(object):
    """
    This is the base class that manages an array of stations that will be
    used during a simulation. It acts like a container that vectorizes the
    station properties to speed up calculations.
    """

    def __init__(self, n):
        self.num_stations = n
        self.x = np.empty(n)
        self.y = np.empty(n)
        self.azimuth = np.empty(n)
        self.elevation = np.empty(n)
        self.height = np.empty(n)
        self.indoor = np.zeros(n, dtype=bool)
        self.active = np.ones(n, dtype=bool)
        self.tx_power = np.empty(n)
        self.rx_power = np.empty(n)
        self.rx_interference = np.empty(n)
        self.ext_interference = np.empty(n)
        self.antenna = np.empty(n, dtype=Antenna)
        self.bandwidth = np.empty(n)
        self.noise_figure = np.empty(n)
        self.noise_temperature = np.empty(n)
        self.thermal_noise = np.empty(n)
        self.total_interference = np.empty(n)
        self.snr = np.empty(n)
        self.sinr = np.empty(n)
        self.sinr_ext = np.empty(n)
        self.inr = np.empty(n)
        self.pfd = np.empty(n)
        self.station_type = StationType.NONE

    def get_station_list(self, id=None) -> list:
        if(id is None):
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
        distance = np.empty([self.num_stations, station.num_stations])
        for i in range(self.num_stations):
            distance[i] = np.sqrt(np.power(self.x[i] - station.x, 2) +
                           np.power(self.y[i] - station.y, 2))
        return distance

    def get_3d_distance_to(self, station) -> np.array:
        distance = np.empty([self.num_stations, station.num_stations])
        for i in range(self.num_stations):
            distance[i] = np.sqrt(np.power(self.x[i] - station.x, 2) +
                           np.power(self.y[i] - station.y, 2) +
                            np.power(self.height[i] - station.height, 2))
        return distance

    def get_elevation_angle(self, station, sat_params) -> dict:
        free_space_angle = np.empty(self.num_stations)
        angle = np.empty(self.num_stations)
        for i in range(self.num_stations):
            # calculate free-space elevation angle according to Attachment A
            rel_x = station.x - self.x[i]
            rel_y = station.y - self.y[i]
            rel_z = station.height - self.height[i]

            gts = np.sqrt(rel_x**2 + rel_y**2)
            theta_0 = np.arctan2(rel_z, gts) # free-space elevation angle
            free_space_angle[i] = np.degrees(theta_0)

            ##
            # calculate apparent elevation angle according to ITU-R P619, Attachment B

            tau_fs1 = 1.728 + 0.5411 * theta_0 + 0.03723 * theta_0**2
            tau_fs2 = 0.1815 + 0.06272 * theta_0 + 0.01380 * theta_0**2
            tau_fs3 = 0.01727 + 0.008288 * theta_0

            # change in elevation angle due to refraction
            tau_fs_deg = 1/(tau_fs1 + sat_params.altitude*tau_fs2 +
                            sat_params.altitude**2*tau_fs3)
            tau_fs = tau_fs_deg / 180. * np.pi

            angle[i] = np.degrees(theta_0 + tau_fs)

        return{'free_space': free_space_angle, 'apparent': angle}

    def get_pointing_vector_to(self, station) -> tuple:

        point_vec_x = station.x- self.x[:,np.newaxis]
        point_vec_y = station.y - self.y[:,np.newaxis]
        point_vec_z = station.height - self.height[:,np.newaxis]

        dist = self.get_3d_distance_to(station)

        phi = np.array(np.rad2deg(np.arctan2(point_vec_y,point_vec_x)),ndmin=2)
        theta = np.rad2deg(np.arccos(point_vec_z/dist))

        return phi, theta

    def get_off_axis_angle(self, station) -> np.array:
        """
        Calculates the off-axis angle between this station and the input station
        """
        Az, b = self.get_pointing_vector_to(station)
        Az0 = self.azimuth

        a = 90 - self.elevation
        C = Az0 - Az
        
        phi = np.arccos(np.cos(np.radians(a))*np.cos(np.radians(b)) \
                        + np.sin(np.radians(a))*np.sin(np.radians(b))*np.cos(np.radians(C)))
        phi_deg = np.degrees(phi)

        return phi_deg
        