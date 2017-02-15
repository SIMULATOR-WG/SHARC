# -*- coding: utf-8 -*-
"""
Created on Wed Feb  8 15:51:29 2017

@author: edgar
"""
import numpy as np

from station_manager import StationManager
from antenna import Antenna
from parameters.parameters_imt import ParametersImt
from macrocell_topology import MacrocellTopology

class DownlinkManager(object):
    
    def __init__(self, k, n_cluster=1):
        self.reset(n_cluster, k)
        
    def reset(self, k, n_cluster=1):
        n_bs = 19
        self.__transmitter = StationManager(n_cluster*n_bs)
        self.__receiver = StationManager(n_cluster*n_bs*k*10)
        self.__link_tx = np.empty(n_cluster*n_bs*k)
        self.__link_rx = np.empty(n_cluster*n_bs*k)
        self.__coupling_loss = np.empty(n_cluster*n_bs*k)
        self.__acir = np.empty(n_cluster*n_bs*k)
        self.__thermal_noise = np.empty(n_cluster*n_bs*k)
        self.__bandwidth = 0
        self.__frequency = 0
        
        self.__topology = MacrocellTopology(ParametersImt.intersite_distance,
                                            ParametersImt.num_clusters)
        
    def generate_base_stations(self, intersite_distance, num_clusters=1):
        (x_coord, y_coord) = self.__topology.get_coordinates()
        # now we set the coordinates
        self.__transmitter.set_x(x_coord)
        self.__transmitter.set_y(y_coord)
        num_stations = self.__transmitter.get_num_stations()
        self.__transmitter.set_height(30*np.ones(num_stations))
        self.__transmitter.set_tx_power(30*np.ones(num_stations))
        self.__transmitter.set_tx_antenna(
            np.array([Antenna(20) for i in range(num_stations)]))