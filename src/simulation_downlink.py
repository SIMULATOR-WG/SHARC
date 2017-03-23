# -*- coding: utf-8 -*-
"""
Created on Wed Jan 11 19:06:41 2017

@author: edgar
"""

import numpy as np
import math

from simulation import Simulation
from parameters.parameters_imt import ParametersImt
from station_manager import StationManager
from topology_macrocell import TopologyMacrocell
from topology_single_base_station import TopologySingleBaseStation
from propagation_free_space import PropagationFreeSpace
from results import Results

class SimulationDownlink(Simulation):
    """
    Implements the flowchart of simulation downlink method
    """
    
    def __init__(self, param: ParametersImt):
        super(SimulationDownlink, self).__init__()
        np.random.seed(0)
        
        self.__param = param
        
        #self.__topology = TopologyMacrocell(self.__param.intersite_distance, self.__param.num_clusters)
        self.__topology = TopologySingleBaseStation(self.__param.intersite_distance/2, self.__param.num_clusters)
        self.__propagation = PropagationFreeSpace()
        
        self.__num_transmitters = self.__param.num_clusters*self.__param.num_base_stations
        self.__num_receivers = self.__param.num_clusters*self.__param.num_base_stations \
                                 *self.__param.ue_k*self.__param.ue_k_m
        
        self.__coupling_loss = np.empty([self.__num_transmitters, self.__num_receivers])
        self.__transmitter = np.empty(self.__num_transmitters)
        self.__receiver = np.empty(self.__num_receivers)            
        self.__link = dict([(bs,list()) for bs in range(self.__num_transmitters)])
        self.__results = Results()
        
    def initialize(self, *args, **kwargs):
        pass
        self.__transmitter = \
            StationManager.generate_imt_base_stations(self.__param,
                                                      self.__topology)
    
    def snapshot(self, *args, **kwargs):
        if not self.__param.static_base_stations:
            # TODO: include support for randomly located base stations by
            # creating the RandomTopology(Topology) class
            pass
        if self.__param.interfered_with:
            self._create_ue()
            self._calculate_coupling_loss()
            self._connect_ue_to_bs()
            self._select_ue()
            self._scheduler()
            self._apply_power_control()
            self._calculate_sinr()
            #self._add_external_interference()
            #self._recalculate_sinr()
            #self._calculate_imt_degradation()
        else:
            #self._beamforming()
            #self._apply_power_control()
            #self._scheduler()
            #self._select_active_bs()
            #self._calculate_other_interference()
            #self._calculate_other_degradation()
            pass
        self._collect_results()

    def finalize(self, *args, **kwargs):
        self.notify_observers(source=__name__, results=self.__results)
        
    def _create_ue(self):
        self.__receiver = StationManager.generate_imt_ue(self.__param,
                                                         self.__topology)    
    
    def _calculate_coupling_loss(self):
        # Calculate distance from transmitters to receivers. The result is a
        # num_tx x num_rx array 
        d = self.__transmitter.get_distance_to(self.__receiver)
        path_loss = self.__propagation.get_loss(distance=d, frequency=self.__param.frequency)
        tx_antenna = self.__transmitter.tx_antenna.astype('float')
        rx_antenna = self.__receiver.rx_antenna.astype('float')
        # replicate columns to have the appropriate size
        tx_gain = np.transpose(np.tile(tx_antenna, (self.__receiver.num_stations, 1)))
        rx_gain = np.tile(rx_antenna, (self.__transmitter.num_stations, 1))
        # calculate coupling loss
        self.__coupling_loss = path_loss - tx_gain - rx_gain
#        self.__coupling_loss = np.maximum(path_loss - tx_gain - rx_gain, 
#                                          ParametersImt.mcl)
        
    def _connect_ue_to_bs(self):
        """
        Link the UE randomly to a BS to which the path coupling loss is within 
        the minimum coupling loss plus the HO margin.
        """
        for ue in range(self.__num_receivers):
            # for each receiver, select the base stations 
            bs_all = np.where(self.__coupling_loss[:,ue] < self.__param.mcl + self.__param.ho_margin)[0]
            bs = np.random.choice(bs_all)
            self.__link[bs].append(ue)
            
    def _select_ue(self):
        """
        Select K UEs randomly from all the UEs linked to one BS as “chosen” 
        UEs. These K “chosen” UEs will be scheduled during this snapshot.
        """
        # TO BE IMPLEMENTED
        pass
        
    def _scheduler(self):
        for bs in range(self.__num_transmitters):
            ue = self.link[bs]
            bandwidth = self.__param.bandwidth / len(ue)
            self.receiver.bandwidth[ue] = bandwidth
    
    def _apply_power_control(self):
        """
        Apply downlink power control algorithm
        """
        # Currently, the maximum transmit power of the base station is equaly
        # divided among the selected UEs
        p_max = math.pow(10, 0.1*self.__param.bs_tx_power)
        # calculate tansmit powers to have a structure such as 
        # {bs_1: [pwr_1, pwr_2,...], ...}, where bs_1 is the base station id, 
        # pwr_1 is the transmit power from bs_1 to ue_1, pwr_2 is the transmit
        # power from bs_1 to ue_2, etc
        self.transmitter.tx_power = dict([(bs,[10*math.log10(p_max/len(self.link[bs]))] * len(self.link[bs])) for bs in range(self.__num_transmitters)])
        
    def _calculate_sinr(self):
        bs_all = [b for b in range(self.__num_transmitters)]
        
        # calculate received power for each UE
        for bs, ue_list in self.link.items():
            self.receiver.rx_power[ue_list] = np.array(self.transmitter.tx_power[bs]) \
                                    - np.array(self.coupling_loss[bs,ue_list])
            # create a list with base stations that generate interference in ue_list
            bs_interf = [b for b in bs_all if b not in [bs]]

            # calculate intra system interference
            for ue in ue_list:
                for bi in bs_interf:
                    interference = self.transmitter.tx_power[bi] \
                        - self.coupling_loss[bi,ue]
                    self.receiver.rx_interference[ue] = 10*math.log10( \
                          math.pow(10, 0.1*self.receiver.rx_interference[ue]) \
                        + np.sum(np.power(10, 0.1*interference)))
        
        self.receiver.thermal_noise = \
            10*math.log10(self.__param.BOLTZMANN_CONSTANT*self.__param.noise_temperature) + \
            10*np.log10(self.receiver.bandwidth * 1e6) + \
            self.receiver.noise_figure
            
        self.receiver.total_interference = \
            10*np.log10(np.power(10, 0.1*self.receiver.rx_interference) + \
                        np.power(10, 0.1*self.receiver.thermal_noise))
        self.receiver.sinr = self.receiver.rx_power - self.receiver.total_interference  
        self.receiver.snr = self.receiver.rx_power - self.receiver.thermal_noise  
                    
    def _collect_results(self):
        self.__results.add_coupling_loss_dl( \
            np.reshape(self.coupling_loss, self.__num_transmitters*self.__num_receivers).tolist())
        
        for bs in range(self.__num_transmitters):
            self.__results.add_tx_power_dl(self.transmitter.tx_power[bs])
            
        self.__results.add_sinr_dl(self.receiver.sinr.tolist())
        self.__results.add_snr_dl(self.receiver.snr.tolist())
        
    @property
    def transmitter(self):
        return self.__transmitter
        
    @property
    def receiver(self):
        return self.__receiver
        
    @property
    def coupling_loss(self):
        return self.__coupling_loss
        
    @property
    def link(self):
        return self.__link
    