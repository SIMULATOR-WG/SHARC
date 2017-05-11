# -*- coding: utf-8 -*-
"""
Created on Fri Feb  3 15:29:48 2017

@author: edgar
"""

import numpy as np

from sharc.station import Station
from sharc.antenna.antenna import Antenna

class StationManager(object):
    """
    This is the base class that manages an array of stations that will be
    used during a simulation. It acts like a container that vectorizes the
    station properties to speed up calculations.
    """

    def __init__(self, n):
        self.__num_stations = n
        self.__x = np.empty(n)
        self.__y = np.empty(n)
        self.__height = np.empty(n)
        self.__active = np.ones(n, dtype=bool)
        self.__tx_power = np.empty(n)
        self.__rx_power = np.empty(n)
        self.__rx_interference = np.empty(n)
        self.__tx_antenna = np.array([Antenna() for i in range(n)])
        self.__rx_antenna = np.array([Antenna() for i in range(n)])
        self.__bandwidth = np.empty(n)
        self.__noise_figure = np.empty(n)
        self.__noise_temperature = np.empty(n)
        self.__thermal_noise = np.empty(n)
        self.__total_interference = np.empty(n)
        self.__snr = np.empty(n)
        self.__sinr = np.empty(n)
        self.__inr = np.empty(n)

    def get_station_list(self,id=None) -> list:
        if(id is None):
            id = range(self.num_stations)
        station_list = list()
        for i in id:
            station_list.append(self.get_station(i))
        return station_list

    def get_station(self,id) -> Station:
        station = Station()
        station.id = id
        station.x = self.x[id]
        station.y = self.y[id]
        station.height = self.height[id]
        station.active = self.active[id]
        station.tx_power = self.tx_power[id]
        station.rx_power = self.rx_power[id]
        station.tx_antenna = self.tx_antenna[id]
        station.rx_antenna = self.rx_antenna[id]
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
            free_space_angle[i] = theta_0

            ##
            # calculate apparent elevation angle according to Attachment B

            tau_fs1 = 1.728 + 0.5411 * theta_0 + 0.03723 * theta_0**2
            tau_fs2 = 0.1815 + 0.06272 * theta_0 + 0.01380 * theta_0**2
            tau_fs3 = 0.01727 + 0.008288 * theta_0

            # change in elevation angle due to refraction
            tau_fs_deg = 1/(tau_fs1 + sat_params.sat_altitude*tau_fs2 +
                            sat_params.sat_altitude**2*tau_fs3)
            tau_fs = tau_fs_deg / 180. * np.pi

            angle[i] = theta_0 + tau_fs

        return{'free_space': free_space_angle, 'apparent': angle}


    @property
    def num_stations(self):
        return self.__num_stations

    @num_stations.setter
    def num_stations(self, value):
        self.__num_stations = value

    @property
    def x(self):
        return self.__x

    @x.setter
    def x(self, value):
        self.__x = np.array(value)

    @property
    def y(self):
        return self.__y

    @y.setter
    def y(self, value):
        self.__y = np.array(value)

    @property
    def height(self):
        return self.__height

    @height.setter
    def height(self, value):
        self.__height = np.array(value)

    @property
    def active(self):
        return self.__active

    @active.setter
    def active(self, value):
        self.__active = np.array(value)

    @property
    def tx_power(self):
        return self.__tx_power

    @tx_power.setter
    def tx_power(self, value):
        self.__tx_power = value

    @property
    def rx_power(self):
        return self.__rx_power

    @rx_power.setter
    def rx_power(self, value):
        self.__rx_power = np.array(value)

    @property
    def rx_interference(self):
        return self.__rx_interference

    @rx_interference.setter
    def rx_interference(self, value):
        self.__rx_interference = np.array(value)

    @property
    def tx_antenna(self):
        return self.__tx_antenna

    @tx_antenna.setter
    def tx_antenna(self, value):
        self.__tx_antenna = np.array(value)

    @property
    def rx_antenna(self):
        return self.__rx_antenna

    @rx_antenna.setter
    def rx_antenna(self, value):
        self.__rx_antenna = np.array(value)

    @property
    def bandwidth(self):
        return self.__bandwidth

    @bandwidth.setter
    def bandwidth(self, value):
        self.__bandwidth = np.array(value)

    @property
    def noise_figure(self):
        return self.__noise_figure

    @noise_figure.setter
    def noise_figure(self, value):
        self.__noise_figure = np.array(value)

    @property
    def noise_temperature(self):
        return self.__noise_temperature

    @noise_temperature.setter
    def noise_temperature(self, value):
        self.__noise_temperature = np.array(value)

    @property
    def thermal_noise(self):
        return self.__thermal_noise

    @thermal_noise.setter
    def thermal_noise(self, value):
        self.__thermal_noise = np.array(value)

    @property
    def total_interference(self):
        return self.__total_interference

    @total_interference.setter
    def total_interference(self, value):
        self.__total_interference = np.array(value)

    @property
    def sinr(self):
        return self.__sinr

    @sinr.setter
    def sinr(self, value):
        self.__sinr = np.array(value)

    @property
    def snr(self):
        return self.__snr

    @snr.setter
    def snr(self, value):
        self.__snr = np.array(value)

    @property
    def inr(self):
        return self.__inr

    @inr.setter
    def inr(self, value):
        self.__inr = np.array(value)
