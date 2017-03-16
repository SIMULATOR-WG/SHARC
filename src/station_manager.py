# -*- coding: utf-8 -*-
"""
Created on Fri Feb  3 15:29:48 2017

@author: edgar
"""

import numpy as np

from station import Station
from antenna import Antenna
from parameters.parameters_imt import ParametersImt

class StationManager(object):
    """
    This is the base class that manages an array of stations that will be 
    used during a simulation. It acts like a container that vectorizes the 
    station properties to speed up calculations. 
    
    TODO: extract the static methods to a factory class
    """
    
    def __init__(self, n):
        self.__num_stations = n
        self.__x = np.empty(n)
        self.__y = np.empty(n)
        self.__height = np.empty(n)
        self.__active = np.ones(n, dtype=bool)
        self.__tx_power = np.empty(n)
        self.__rx_power = np.empty(n)
        self.__tx_antenna = np.array([Antenna() for i in range(n)])
        self.__rx_antenna = np.array([Antenna() for i in range(n)])
        
    @staticmethod
    def generate_imt_base_stations(num_base_stations, topology):
        imt_base_stations = StationManager(num_base_stations)
        (x_coord, y_coord) = topology.get_coordinates()
        # now we set the coordinates
        imt_base_stations.set_x(x_coord)
        imt_base_stations.set_y(y_coord)
        imt_base_stations.set_height(ParametersImt.bs_height*np.ones(num_base_stations))
        imt_base_stations.set_active(np.ones(num_base_stations, dtype=bool))
        imt_base_stations.set_tx_power(ParametersImt.bs_tx_power*np.ones(num_base_stations))
        imt_base_stations.set_tx_antenna(
            np.array([Antenna(ParametersImt.bs_tx_antenna_gain) for i in range(num_base_stations)]))
        imt_base_stations.set_rx_antenna(
            np.array([Antenna(ParametersImt.bs_rx_antenna_gain) for i in range(num_base_stations)]))                
        return imt_base_stations
        
    @staticmethod
    def generate_imt_ue(num_ue, x_limits, y_limits):
        imt_ue = StationManager(num_ue)
        x_coord = (x_limits[1]-x_limits[0])*np.random.random(num_ue)+x_limits[0]
        y_coord = (y_limits[1]-y_limits[0])*np.random.random(num_ue)+y_limits[0]
        imt_ue.set_x(x_coord)
        imt_ue.set_y(y_coord)
#        imt_ue.set_x(300)
#        imt_ue.set_y(400)
        imt_ue.set_height(ParametersImt.ue_height*np.ones(num_ue))
        imt_ue.set_tx_power(ParametersImt.ue_tx_power*np.ones(num_ue))
        imt_ue.set_tx_antenna(
            np.array([Antenna(ParametersImt.ue_tx_antenna_gain) for i in range(num_ue)]))
        imt_ue.set_rx_antenna(
            np.array([Antenna(ParametersImt.ue_rx_antenna_gain) for i in range(num_ue)]))                
        return imt_ue
        
    def get_station_list(self,id=None):
        if(id is None):
            id = range(self.num_stations)
        station_list = list()
        for i in id:
            station_list.append(self.get_station(i))
        return station_list
      
    def get_station(self,id):
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
        
    @property
    def num_stations(self):
        return self.__num_stations
        
    @num_stations.setter
    def num_stations(self, value):
        self.__num_stations = np.array(value)
        
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
        self.__tx_power = np.array(value)
    
    @property
    def rx_power(self):
        return self.__rx_power
        
    @rx_power.setter
    def rx_power(self, value):
        self.__rx_power = np.array(value)
            
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
        
