# -*- coding: utf-8 -*-
"""
Created on Wed Jan 11 19:06:41 2017

@author: edgar
"""

import numpy as np

from simulation import Simulation
from parameters.parameters_imt import ParametersImt
from station_manager import StationManager
from macrocell_topology import MacrocellTopology

class SimulationDownlink(Simulation):
    """
    Implements the flowchart of simulation downlink method
    """
    
    def __init__(self):
        self.__bs = ParametersImt.num_base_stations
        self.__c = ParametersImt.num_clusters
        self.__k = ParametersImt.ue_k
        self.__m = ParametersImt.ue_k_m
        
        self.__topology = MacrocellTopology(ParametersImt.intersite_distance,
                                    ParametersImt.num_clusters)

        self.reset()
    
    def reset(self):
        #self.__receiver = StationManager(c*bs*k*m)
        self.__link_tx = np.empty(self.__c*self.__bs*self.__k)
        self.__link_rx = np.empty(self.__c*self.__bs*self.__k)
        self.__coupling_loss = np.empty(self.__c*self.__bs*self.__k)
        self.__acir = np.empty(self.__c*self.__bs*self.__k)
        self.__thermal_noise = np.empty(self.__c*self.__bs*self.__k)
        self.__bandwidth = ParametersImt.bandwidth
        self.__frequency = ParametersImt.frequency        
        
    def initialize(self, *args, **kwargs):
        self.__transmitter = \
            StationManager.generate_imt_base_stations(self.__c*self.__bs, self.__topology)
    
    def snapshot(self, *args, **kwargs):
        if not ParametersImt.static_base_stations:
            # TODO: include support for randomly located base stations by
            # creating the RandomTopology(Topology) class
            pass
        if ParametersImt.interfered_with:
            #self.create_ue()
            self.__receiver = StationManager.generate_imt_ue(self.c*self.bs*self.__k*self.__m,
                                                             self.__topology.get_x_limits(),
                                                             self.__topology.get_y_limits())
            self.calculate_coupling_loss()
            self.connect_ue_to_bs()
            #self.apply_power_control()
            #self.scheduler()
            self.calculate_powers()
            self.add_external_interference()
            self.recalculate_sinr()
            self.calculate_imt_degradation()
        else:
            #self.beamforming()
            #self.apply_power_control()
            #self.scheduler()
            #self.select_active_bs()
            #self.calculate_other_interference()
            #self.calculate_other_degradation()
            pass
        #self.collect_results()

    def finalize(self, *args, **kwargs):
        pass
        
