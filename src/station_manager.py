# -*- coding: utf-8 -*-
"""
Created on Fri Feb  3 15:29:48 2017

@author: edgar
"""

import numpy as np

from station import Station
from antenna import Antenna
from parameters.ParametersImt import ParametersImt

class StationManager(object):
    """
    This is the base class that manages an array of stations that will be 
    used during a simulation. It acts like a container that vectorizes the 
    station properties to speed up calculations. 
    
    TODO: extract the static methods to a factory class
    """
    
    def __init__(self, n):
        self.reset(n)
        
    def reset(self, n):
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
        x_coord = (x_limits(1)-x_limits(0))*np.random.random(num_ue)+x_limits(0)
        y_coord = (y_limits(1)-y_limits(0))*np.random.random(num_ue)+y_limits(0)
        imt_ue.set_x(x_coord)
        imt_ue.set_x(y_coord)
        imt_ue.set_height(ParametersImt.ue_height*np.ones(num_ue))
        imt_ue.set_tx_power(ParametersImt.ue_tx_power*np.ones(num_ue))
        imt_ue.set_tx_antenna(
            np.array([Antenna(ParametersImt.ue_tx_antenna_gain) for i in range(num_ue)]))
        imt_ue.set_rx_antenna(
            np.array([Antenna(ParametersImt.ue_rx_antenna_gain) for i in range(num_ue)]))                
        return imt_ue
        
    def get_station_list(self,id=None):
        if(id is None):
            id = range(self.get_num_stations())
        station_list = list()
        for i in id:
            station_list.append(self.get_station(i))
        return station_list
      
    def get_station(self,id):
        station = Station()
        station.set_id(id)
        station.set_x(self.get_x(id))
        station.set_y(self.get_y(id))
        station.set_height(self.get_height(id))
        station.set_active(self.get_active(id))
        station.set_tx_power(self.get_tx_power(id))
        station.set_rx_power(self.get_rx_power(id))
        station.set_tx_antenna(self.get_tx_antenna(id))
        station.set_rx_antenna(self.get_rx_antenna(id))
        return station
        
    def get_num_stations(self):
        return self.__num_stations
        
    def set_x(self, x, id=None):
        if(id is None):
            id = range(self.get_num_stations())
        self.__x[id] = np.asarray(x)

    def get_x(self, id=None):
        if(id is None):
            id = range(self.get_num_stations())
        return self.__x[id]

    def set_y(self, y, id=None):
        if(id is None):
            id = range(self.get_num_stations())
        self.__y[id] = np.asarray(y)

    def get_y(self, id=None):
        if(id is None):
            id = range(self.get_num_stations())
        return self.__y[id]

    def set_height(self, height, id=None):
        if(id is None):
            id = range(self.get_num_stations())
        self.__height[id] = np.asarray(height)

    def get_height(self, id=None):
        if(id is None):
            id = range(self.get_num_stations())
        return self.__height[id]

    def set_active(self, active, id=None):
        if(id is None):
            id = range(self.get_num_stations())
        self.__active[id] = np.asarray(active)

    def get_active(self, id=None):
        if(id is None):
            id = range(self.get_num_stations())
        return self.__active[id]

    def set_tx_power(self, tx_power, id=None):
        if(id is None):
            id = range(self.get_num_stations())
        self.__tx_power[id] = np.asarray(tx_power)

    def get_tx_power(self, id=None):
        if(id is None):
            id = range(self.get_num_stations())
        return self.__tx_power[id]

    def set_rx_power(self, rx_power, id=None):
        if(id is None):
            id = range(self.get_num_stations())
        self.__rx_power[id] = np.asarray(rx_power)

    def get_rx_power(self, id=None):
        if(id is None):
            id = range(self.get_num_stations())
        return self.__rx_power[id]

    def set_tx_antenna(self, tx_antenna, id=None):
        if(id is None):
            id = range(self.get_num_stations())
        self.__tx_antenna[id] = tx_antenna

    def get_tx_antenna(self, id=None):
        if(id is None):
            id = range(self.get_num_stations())
        return self.__tx_antenna[id]

    def set_rx_antenna(self, rx_antenna, id=None):
        if(id is None):
            id = range(self.get_num_stations())
        self.__rx_antenna[id] = rx_antenna

    def get_rx_antenna(self, id=None):
        if(id is None):
            id = range(self.get_num_stations())
        return self.__rx_antenna[id]
