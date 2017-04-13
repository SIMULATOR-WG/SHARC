# -*- coding: utf-8 -*-
"""
Created on Fri Apr  7 17:02:35 2017

@author: edgar
"""

import numpy as np
import random
import math

from sharc.simulation import Simulation
from sharc.parameters.parameters_imt import ParametersImt
from sharc.station_factory import StationFactory
from sharc.topology.topology_factory import TopologyFactory
from sharc.propagation.propagation_free_space import PropagationFreeSpace
from sharc.results import Results

class SimulationUplink(Simulation):
    """
    Implements the flowchart of simulation downlink method
    """
    
    def __init__(self, param: ParametersImt):
        super(SimulationUplink, self).__init__()
        self.param = param
        self.topology = TopologyFactory.createTopology(self.param)
        self.propagation = PropagationFreeSpace()
        
        num_ue = self.param.num_clusters*self.param.num_base_stations \
                                 *self.param.ue_k*self.param.ue_k_m
        num_bs = self.param.num_clusters*self.param.num_base_stations
        
        self.coupling_loss = np.empty([num_bs, num_ue])
        
        self.ue = np.empty(num_ue)
        self.bs = np.empty(num_bs)
        self.other = np.empty(1)
        
        # this attribute indicates the list of UE's that are connected to each
        # base station. The position the the list indicates the resource block
        # group that is allocated to the given UE
        self.link = dict([(bs,list()) for bs in range(num_bs)])
        
        # calculates the number of RB per BS
        self.num_rb_per_bs = math.trunc((1-self.param.guard_band_ratio)* \
                            self.param.bandwidth /self.param.rb_bandwidth)
        # calculates the number of RB per UE on a given BS
        self.num_rb_per_ue = math.trunc(self.num_rb_per_bs/self.param.ue_k)
        
        self.results = Results()
        
    def initialize(self, *args, **kwargs):
        self.bs = StationFactory.generate_imt_base_stations(self.param,
                                                            self.topology)
    
    def snapshot(self, *args, **kwargs):
        self.create_ue()
        self.calculate_coupling_loss()
        self.connect_ue_to_bs()
        self.select_ue()
        self.scheduler()
        self.power_control()        
        if not self.param.static_base_stations:
            # TODO: include support for randomly located base stations by
            # creating the RandomTopology(Topology) class
            pass
        if self.param.interfered_with:
            pass
            #self.calculate_sinr()
            #self.add_external_interference()
            #self.recalculate_sinr()
            #self.calculate_imt_degradation()
        else:
            self.calculate_external_interference()
            #self.calculate_external_degradation()
            pass
        self.collect_results()

    def finalize(self, *args, **kwargs):
        self.notify_observers(source=__name__, results=self.results)
        
    def create_ue(self):
        self.ue = StationFactory.generate_imt_ue(self.param, self.topology)    
    
    def calculate_coupling_loss(self):
        """
        Calculates the path coupling loss from each UE to all BS's. Result is 
        stored as a numpy array (self.coupling_loss) with dimensions 
        num_bs x num_ue 
        """
        # Calculate distance from transmitters to receivers. The result is a
        # num_bs x num_ue array 
        d = self.bs.get_distance_to(self.ue)
        path_loss = self.propagation.get_loss(distance=d, frequency=self.param.frequency)
        bs_antenna = self.bs.rx_antenna.astype('float')
        ue_antenna = self.ue.tx_antenna.astype('float')
        # replicate columns to have the appropriate size
        bs_gain = np.transpose(np.tile(bs_antenna, (self.ue.num_stations, 1)))
        ue_gain = np.tile(ue_antenna, (self.bs.num_stations, 1))
        # calculate coupling loss
        self.coupling_loss = path_loss - bs_gain - ue_gain
#        self.coupling_loss = np.maximum(path_loss - tx_gain - rx_gain, 
#                                          ParametersImt.mcl)
        
    def connect_ue_to_bs(self):
        """
        Link the UE randomly to a BS to which the path coupling loss is within 
        the minimum coupling loss plus the HO margin.
        """
        for ue in range(self.ue.num_stations):
            # for each receiver, select the base stations 
            bs_all = np.where(self.coupling_loss[:,ue] < self.param.mcl + self.param.ho_margin)[0]
            bs = np.random.choice(bs_all)
            self.link[bs].append(ue)
            
    def select_ue(self):
        """
        Select K UEs randomly from all the UEs linked to one BS as “chosen” 
        UEs. These K “chosen” UEs will be scheduled during this snapshot.
        """
        for bs in range(self.bs.num_stations):
            # select K UE's among the ones that are connected to BS
            random.shuffle(self.link[bs])
            K = self.param.ue_k
            #K = np.random.choice(range(1, len(self.link[bs])))
            del self.link[bs][K:]
            # Activate the selected UE's
            self.ue.active[self.link[bs]] = np.ones(K, dtype=bool)
        
    def scheduler(self):
        """
        This scheduler divides the available resource blocks among UE's for 
        a given BS
        """
        for bs in range(self.bs.num_stations):
            ue_list = self.link[bs]
            self.ue.bandwidth[ue_list] = self.num_rb_per_ue*self.param.rb_bandwidth
    
    def power_control(self):
        """
        Apply uplink power control algorithm
        """
        # Currently, there is no power control; then, set it to maximum
        self.ue.tx_power = self.param.ue_tx_power*np.ones(self.ue.num_stations)
        
    def calculate_sinr(self):
        """
        Calculates the uplink SINR for each UE. This is useful only in the 
        cases when IMT system is interfered by other system
        TODO: not working yet
        """
        pass
#        bs_all = [b for b in range(self.bs.num_stations)]
#        self.bs.rx_power = dict([(bs,[-300] * len(self.link[bs])) for bs in range(self.bs.num_stations)])
#        
#        # calculate uplink received power for each BS
#        for bs, ue_list in self.link.items():
#            self.bs.rx_power[bs] = np.array(self.ue.tx_power[ue_list]) \
#                                    - np.array(self.coupling_loss[bs,ue_list])
#            # create a list with base stations that generate interference in ue_list
#            bs_interf = [b for b in bs_all if b not in [bs]]
#
#            # calculate intra system interference
#            for bi in bs_interf:
#                ui = self.link[bi]
#                interference = self.ue.tx_power[ui] - self.coupling_loss[bi,ui]
#                self.bs.rx_interference[ue] = 10*math.log10( \
#                    math.pow(10, 0.1*self.bs.rx_interference[ue]) 
#                    + np.sum(np.power(10, 0.1*interference)))
#        
#        self.bs.thermal_noise = \
#            10*math.log10(self.param.BOLTZMANN_CONSTANT*self.param.noise_temperature) + \
#            10*np.log10(self.bs.bandwidth * 1e6) + \
#            self.bs.noise_figure
#            
#        self.bs.total_interference = \
#            10*np.log10(np.power(10, 0.1*self.bs.rx_interference) + \
#                        np.power(10, 0.1*self.bs.thermal_noise))
#        self.bs.sinr = self.bs.rx_power - self.bs.total_interference  
#        self.bs.snr = self.bs.rx_power - self.bs.thermal_noise  
                    
    def calculate_external_interference(self):
        pass


    def collect_results(self):
        self.results.add_coupling_loss_dl( \
            np.reshape(self.coupling_loss, self.num_ue*self.num_bs).tolist())

        for bs in range(self.bs.num_stations):
            self.results.add_tx_power_dl(self.ue.tx_power[bs])
            
        # select the active stations
        ids = np.where(self.bs.active)

        self.results.add_sinr_dl(self.bs.sinr[ids].tolist())
        self.results.add_snr_dl(self.bs.snr[ids].tolist())
        
    