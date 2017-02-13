# -*- coding: utf-8 -*-
"""
Created on Wed Feb  8 15:51:29 2017

@author: edgar
"""
import numpy as np
import math

from station_manager import StationManager
from antenna import Antenna

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
        
    def generate_base_stations(self, intersite_distance, num_clusters=1):
        d = intersite_distance
        h = (d/3)*math.sqrt(3)/2
        # these are the coordinates of the central cluster
        x_central = np.array([0, d, d/2, -d/2, -d, -d/2, 
                         d/2, 2*d, 3*d/2, d, 0, -d, 
                         -3*d/2, -2*d, -3*d/2, -d, 0, d, 3*d/2])
        y_central = np.array([0, 0, 3*h, 3*h, 0, -3*h, 
                         -3*h, 0, 3*h, 6*h, 6*h, 6*h, 
                         3*h, 0, -3*h, -6*h, -6*h, -6*h, -3*h])
        x_coord = np.copy(x_central)
        y_coord = np.copy(y_central)
        # other clusters are calculated by shifting the central cluster
        if num_clusters == 7:
            x_shift = np.array([7*d/2, -d/2, -4*d, -7*d/2, d/2, 4*d])
            y_shift = np.array([9*h, 15*h, 6*h, -9*h, -15*h, -6*h])
            for x, y in zip(x_shift, y_shift):
                x_coord = np.concatenate((x_coord, x_central+x))
                y_coord = np.concatenate((y_coord, y_central+y))
        # now we set the coordinates
        self.__transmitter.set_x(x_coord)
        self.__transmitter.set_y(y_coord)
        num_stations = self.__transmitter.get_num_stations()
        self.__transmitter.set_height(30*np.ones(num_stations))
        self.__transmitter.set_tx_power(30*np.ones(num_stations))
        self.__transmitter.set_tx_antenna(
            np.array([Antenna(20) for i in range(num_stations)]))