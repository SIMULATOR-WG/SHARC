# -*- coding: utf-8 -*-
"""
Created on Wed Jan 11 19:06:41 2017

@author: edgar
"""

import numpy as np

from simulation import Simulation
from parameters.parameters_imt import ParametersImt
from station_manager import StationManager
from topology_macrocell import TopologyMacrocell
from topology_single_base_station import TopologySingleBaseStation
from propagation_free_space import PropagationFreeSpace

class SimulationDownlink(Simulation):
    """
    Implements the flowchart of simulation downlink method
    """
    
    def __init__(self):
        np.random.seed(0)
        
        self.__bs = ParametersImt.num_base_stations
        self.__c = ParametersImt.num_clusters
        self.__k = ParametersImt.ue_k
        self.__m = ParametersImt.ue_k_m
        
        #self.__topology = TopologyMacrocell(ParametersImt.intersite_distance, ParametersImt.num_clusters)
        self.__topology = TopologySingleBaseStation(ParametersImt.intersite_distance/2)
        self.__propagation = PropagationFreeSpace()

        self.__coupling_loss = np.empty([self.__c*self.__bs, self.__k*self.__m])
        self.__transmitter = \
            StationManager.generate_imt_base_stations(self.__c*self.__bs, self.__topology)
        self.__bandwidth = ParametersImt.bandwidth
        self.__frequency = ParametersImt.frequency      
      
        
    def initialize(self, *args, **kwargs):
        pass
#        self.__transmitter = \
#            StationManager.generate_imt_base_stations(self.__c*self.__bs, self.__topology)
    
    def snapshot(self, *args, **kwargs):
        if not ParametersImt.static_base_stations:
            # TODO: include support for randomly located base stations by
            # creating the RandomTopology(Topology) class
            pass
        if ParametersImt.interfered_with:
            #self.create_ue()
            self.__receiver = StationManager.generate_imt_ue(self.__c*self.__bs*self.__k*self.__m,
                                                             self.__topology.get_x_limits(),
                                                             self.__topology.get_y_limits())
            self._calculate_coupling_loss()
            self._connect_ue_to_bs()
            #self._apply_power_control()
            #self._scheduler()
            self._calculate_powers()
            self._add_external_interference()
            self._recalculate_sinr()
            self._calculate_imt_degradation()
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
        
    def _calculate_coupling_loss(self):
        # Calculate distance from transmitters to receivers. The result is a
        # num_tx x num_rx array 
        d = self.__transmitter.get_distance_to(self.__receiver)
        path_loss = self.__propagation.get_loss(distance=d, frequency=self.__frequency)
        tx_antenna = self.__transmitter.get_tx_antenna()
        rx_antenna = self.__receiver.get_rx_antenna()
        # replicate columns to have the appropriate size
        tx_gain = np.tile(tx_antenna, (self.__transmitter.get_num_stations(), self.__k))
        rx_gain = np.tile(rx_antenna, (self.__transmitter.get_num_stations(), 1))
        # calculate coupling loss
        self.__coupling_loss = np.maximum(tx_gain + path_loss + rx_gain, 
                                          ParametersImt.mcl)
        
    def _connect_ue_to_bs(self):
        """
        Link the UE randomly to a BS to which the path coupling loss is within 
        the smallest coupling loss plus the HO margin.
        """
        for i in range(self.__receiver.get_num_stations()):
            # for each receiver, select the base stations 
            bs_all = np.where(self.__coupling_loss[:,i] > ParametersImt.mcl + ParametersImt.ho_margin)[0]
            bs_chosen = np.random.choice(bs_all)