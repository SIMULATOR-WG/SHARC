# -*- coding: utf-8 -*-
"""
Created on Wed Jan 11 19:06:41 2017

@author: edgar
"""

import numpy as np
import random
import math
import sys

from sharc.simulation import Simulation
from sharc.parameters.parameters_imt import ParametersImt
from sharc.parameters.parameters_antenna_imt import ParametersAntennaImt
from sharc.station_factory import StationFactory
from sharc.topology.topology_factory import TopologyFactory
from sharc.propagation.propagation_free_space import PropagationFreeSpace
from sharc.propagation.propagation_close_in import PropagationCloseIn
from sharc.results import Results

class SimulationDownlink(Simulation):
    """
    Implements the flowchart of simulation downlink method
    """

    def __init__(self, param: ParametersImt, param_ant: ParametersAntennaImt):
        super(SimulationDownlink, self).__init__()
        np.random.seed(0)

        self.param = param
        self.param_imt_antenna = param_ant

        self.topology = TopologyFactory.createTopology(self.param)

        if self.param.channel_model == "FSPL":
            self.propagation = PropagationFreeSpace()
        elif self.param.channel_model == "CI":
            self.propagation = PropagationCloseIn( self.param.topology,
                                                   self.param.line_of_sight_prob )
        else:
            sys.stderr.write("error: invalid parameter channel_model\n")
            sys.exit(1)


        self.num_transmitters = self.param.num_clusters*self.param.num_base_stations
        self.num_receivers = self.param.num_clusters*self.param.num_base_stations \
                                 *self.param.ue_k*self.param.ue_k_m

        self.coupling_loss = np.empty([self.num_transmitters, self.num_receivers])
        self.transmitter = np.empty(self.num_transmitters)
        self.receiver = np.empty(self.num_receivers)
        self.link = dict([(bs,list()) for bs in range(self.num_transmitters)])
        self.results = Results()

    def initialize(self, *args, **kwargs):
        pass
        self.transmitter = \
            StationFactory.generate_imt_base_stations(self.param,
                                                      self.param_imt_antenna,
                                                      self.topology)

    def snapshot(self, *args, **kwargs):
        if not self.param.static_base_stations:
            # TODO: include support for randomly located base stations by
            # creating the RandomTopology(Topology) class
            pass
        if self.param.interfered_with:
            self.create_ue()
            self.calculate_coupling_loss()
            self.connect_ue_to_bs()
            self.select_ue()
            self.scheduler()
            self.power_control()
            self.calculate_sinr()
            #self.add_external_interference()
            #self.recalculate_sinr()
            #self.calculate_imt_degradation()
        else:
            #self.beamforming()
            #self.apply_power_control()
            #self.scheduler()
            #self.select_active_bs()
            #self.calculate_other_interference()
            #self.calculate_other_degradation()
            pass
        self.collect_results()

    def finalize(self, *args, **kwargs):
        self.notify_observers(source=__name__, results=self.results)

    def create_ue(self):
        self.receiver = StationFactory.generate_imt_ue(self.param,
                                                       self.param_imt_antenna,
                                                       self.topology)

    def calculate_coupling_loss(self):
        # Calculate distance from transmitters to receivers. The result is a
        # num_tx x num_rx array
        d = self.transmitter.get_distance_to(self.receiver)
        path_loss = self.propagation.get_loss(distance=d, frequency=self.param.frequency)
        tx_antenna = self.transmitter.tx_antenna.astype('float')
        rx_antenna = self.receiver.rx_antenna.astype('float')
        # replicate columns to have the appropriate size
        tx_gain = np.transpose(np.tile(tx_antenna, (self.receiver.num_stations, 1)))
        rx_gain = np.tile(rx_antenna, (self.transmitter.num_stations, 1))
        # calculate coupling loss
        self.coupling_loss = path_loss - tx_gain - rx_gain

    def connect_ue_to_bs(self):
        """
        Link the UE randomly to a BS to which the path coupling loss is within
        the minimum coupling loss plus the HO margin.
        """
        for ue in range(self.num_receivers):
            # for each receiver, select the base stations
            minimum_coupling_loss = np.amin(self.coupling_loss[:,ue])
            bs_all = np.where(self.coupling_loss[:,ue] < minimum_coupling_loss + self.param.ho_margin)[0]
            bs = np.random.choice(bs_all)
            self.link[bs].append(ue)

    def select_ue(self):
        """
        Select K UEs randomly from all the UEs linked to one BS as “chosen”
        UEs. These K “chosen” UEs will be scheduled during this snapshot.
        """
        for bs in range(self.num_transmitters):
            # select K UE's among the ones that are connected to BS
            random.shuffle(self.link[bs])
            K = self.param.ue_k
            #K = np.random.choice(range(1, len(self.link[bs])))
            del self.link[bs][K:]
            # Activate the selected UE's
            self.receiver.active[self.link[bs]] = np.ones(K, dtype=bool)

    def scheduler(self):
        """
        This scheduler divides the available bandwidth among UE's ona given BS
        """
        for bs in range(self.num_transmitters):
            ue = self.link[bs]
            bandwidth = self.param.bandwidth / len(ue)
            self.receiver.bandwidth[ue] = bandwidth

    def power_control(self):
        """
        Apply downlink power control algorithm
        """
        # Currently, the maximum transmit power of the base station is equaly
        # divided among the selected UEs
        p_max = math.pow(10, 0.1*self.param.bs_tx_power)
        # calculate tansmit powers to have a structure such as
        # {bs_1: [pwr_1, pwr_2,...], ...}, where bs_1 is the base station id,
        # pwr_1 is the transmit power from bs_1 to ue_1, pwr_2 is the transmit
        # power from bs_1 to ue_2, etc
        self.transmitter.tx_power = dict([(bs,[10*math.log10(p_max/len(self.link[bs]))] * len(self.link[bs])) for bs in range(self.num_transmitters)])

    def calculate_sinr(self):
        bs_all = [b for b in range(self.num_transmitters)]

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
                        math.pow(10, 0.1*self.receiver.rx_interference[ue])
                        + np.sum(np.power(10, 0.1*interference)))

        self.receiver.thermal_noise = \
            10*math.log10(self.param.BOLTZMANN_CONSTANT*self.param.noise_temperature) + \
            10*np.log10(self.receiver.bandwidth * 1e6) + \
            self.receiver.noise_figure

        self.receiver.total_interference = \
            10*np.log10(np.power(10, 0.1*self.receiver.rx_interference) + \
                        np.power(10, 0.1*self.receiver.thermal_noise))
        self.receiver.sinr = self.receiver.rx_power - self.receiver.total_interference
        self.receiver.snr = self.receiver.rx_power - self.receiver.thermal_noise

    def collect_results(self):
        self.results.add_coupling_loss_dl( \
            np.reshape(self.coupling_loss, self.num_transmitters*self.num_receivers).tolist())

        for bs in range(self.num_transmitters):
            self.results.add_tx_power_dl(self.transmitter.tx_power[bs])

        # select the active stations
        ids = np.where(self.receiver.active)

        self.results.add_sinr_dl(self.receiver.sinr[ids].tolist())
        self.results.add_snr_dl(self.receiver.snr[ids].tolist())

