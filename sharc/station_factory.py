# -*- coding: utf-8 -*-
"""
Created on Thu Mar 23 16:37:32 2017

@author: edgar
"""

import numpy as np

from sharc.parameters.parameters_imt import ParametersImt
from sharc.station_manager import StationManager
from sharc.antenna import Antenna
from sharc.topology import Topology

class StationFactory(object):
    
    @staticmethod
    def generate_imt_base_stations(param: ParametersImt, topology: Topology):
        num_bs = param.num_clusters*param.num_base_stations
        imt_base_stations = StationManager(num_bs)
        # now we set the coordinates
        imt_base_stations.x = topology.x
        imt_base_stations.y = topology.y
        imt_base_stations.height = param.bs_height*np.ones(num_bs)
        imt_base_stations.active = np.ones(num_bs, dtype=bool)
        imt_base_stations.tx_power = param.bs_tx_power*np.ones(num_bs)
        imt_base_stations.tx_antenna = \
            np.array([Antenna(param.bs_tx_antenna_gain) for i in range(num_bs)])
        imt_base_stations.rx_antenna = \
            np.array([Antenna(param.bs_rx_antenna_gain) for i in range(num_bs)])  
        imt_base_stations.bandwidth = param.bandwidth*np.ones(num_bs)
        imt_base_stations.noise_figure = param.bs_noise_figure*np.ones(num_bs)
        return imt_base_stations
        
    @staticmethod
    def generate_imt_ue(param: ParametersImt, topology: Topology):
        num_ue = param.num_clusters*param.num_base_stations*param.ue_k*param.ue_k_m
        imt_ue = StationManager(num_ue)
        #imt_ue.x = (topology.x_max - topology.x_min)*np.random.random(num_ue) + topology.x_min
        #imt_ue.y = (topology.y_max - topology.y_min)*np.random.random(num_ue) + topology.y_min
        ue_x = list()
        ue_y = list()
        for bs in range(topology.x.size):
            x_min = topology.x[bs] - topology.cell_radius
            x_max = topology.x[bs] + topology.cell_radius
            y_min = topology.y[bs] - topology.cell_radius
            y_max = topology.y[bs] + topology.cell_radius
            x = (x_max - x_min)*np.random.random(param.ue_k*param.ue_k_m) + x_min
            y = (y_max - y_min)*np.random.random(param.ue_k*param.ue_k_m) + y_min
            ue_x.extend(x)
            ue_y.extend(y)
        imt_ue.x = np.array(ue_x)
        imt_ue.y = np.array(ue_y)   
        imt_ue.active = np.zeros(num_ue, dtype=bool)
        imt_ue.height = param.ue_height*np.ones(num_ue)
        imt_ue.tx_power = param.ue_tx_power*np.ones(num_ue)
        imt_ue.rx_interference = -300*np.ones(num_ue)
        imt_ue.tx_antenna = \
            np.array([Antenna(param.ue_tx_antenna_gain) for i in range(num_ue)])
        imt_ue.rx_antenna = \
            np.array([Antenna(param.ue_rx_antenna_gain) for i in range(num_ue)])   
        imt_ue.bandwidth = param.bandwidth*np.ones(num_ue)
        imt_ue.noise_figure = param.ue_noise_figure*np.ones(num_ue)
        return imt_ue
        
    