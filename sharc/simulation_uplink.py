# -*- coding: utf-8 -*-
"""
Created on Fri Apr  7 17:02:35 2017

@author: edgar
"""

import numpy as np
import random
import math
import sys

from sharc.simulation import Simulation
from sharc.parameters.parameters_imt import ParametersImt
from sharc.parameters.parameters_fss import ParametersFss
from sharc.station_factory import StationFactory
from sharc.station_manager import StationManager
from sharc.topology.topology_factory import TopologyFactory
from sharc.propagation.propagation import Propagation
from sharc.propagation.propagation_free_space import PropagationFreeSpace
from sharc.propagation.propagation_close_in import PropagationCloseIn
from sharc.propagation.propagation_p619 import PropagationP619
from sharc.propagation.propagation_sat_simple import PropagationSatSimple


from sharc.results import Results

class SimulationUplink(Simulation):
    """
    Implements the flowchart of simulation downlink method
    """

    def __init__(self, param: ParametersImt):
        super(SimulationUplink, self).__init__()
        self.param = param
        self.param_system = ParametersFss()

        self.topology = TopologyFactory.createTopology(self.param)

        if self.param.channel_model == "FSPL":
            self.propagation_imt = PropagationFreeSpace()
        elif self.param.channel_model == "CI":
            self.propagation_imt = PropagationCloseIn(self.param.topology,
                                                      self.param.line_of_sight_prob)
        else:
            sys.stderr.write("error: invalid parameter channel_model" + self.param.channel_model
                             + "in IMT propagation\n")
            sys.exit(1)


        if self.param_system.channel_model == "FSPL":
            self.propagation_system = PropagationFreeSpace()
        elif self.param_system.channel_model == "SatelliteSimple":
            self.propagation_system = PropagationSatSimple(self.param_system.line_of_sight_prob)
        elif self.param_system.channel_model == "P619":
            self.propagation_system = PropagationP619()
        else:
            sys.stderr.write("error: invalid parameter channel_model" + self.param_system.channel_model
                             + "in satellite propagation\n")
            sys.exit(1)

        num_ue = self.param.num_clusters*self.param.num_base_stations \
                                 *self.param.ue_k*self.param.ue_k_m
        num_bs = self.param.num_clusters*self.param.num_base_stations

        self.coupling_loss = np.empty([num_bs, num_ue])
        self.coupling_loss_ue_sat = np.empty(num_ue)
        self.coupling_loss_bs_sat = np.empty(num_bs)

        self.ue = np.empty(num_ue)
        self.bs = np.empty(num_bs)
        self.system = np.empty(1)

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
        self.create_system()
        self.create_ue()
        self.coupling_loss = np.transpose( \
                             self.calculate_coupling_loss(self.ue, self.bs,
                                                          self.propagation_imt))
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
        self.results.write_files()
        self.notify_observers(source=__name__, results=self.results)

    def create_ue(self):
        self.ue = StationFactory.generate_imt_ue(self.param, self.topology)

    def create_system(self):
        self.system = StationFactory.generate_fss_stations(self.param_system)

    def calculate_coupling_loss(self,
                                station_a: StationManager,
                                station_b: StationManager,
                                propagation: Propagation) -> np.array:
        """
        Calculates the path coupling loss from each stationA to all stationB.
        Result is returned as a numpy array with dimensions num_a x num_b
        """
        # Calculate distance from transmitters to receivers. The result is a
        # num_bs x num_ue array
        d = station_a.get_3d_distance_to(station_b)

        if station_b.is_satellite:
            elevation_angles = station_a.get_elevation_angle(station_b, self.param_system)
            path_loss = propagation.get_loss(distance=d, frequency=self.param.frequency,
                                             elevation=elevation_angles, sat_params = self.param_system,
                                             earth_to_space = True)
        else:
            path_loss = propagation.get_loss(distance=d, frequency=self.param.frequency)
        antenna_a = station_a.tx_antenna.astype('float')
        antenna_b = station_b.rx_antenna.astype('float')
        # replicate columns to have the appropriate size
        gain_a = np.transpose(np.tile(antenna_a, (station_b.num_stations, 1)))
        gain_b = np.tile(antenna_b, (station_a.num_stations, 1))
        # calculate coupling loss
        return path_loss - gain_a - gain_b
#        self.coupling_loss = np.maximum(path_loss - tx_gain - rx_gain,
#                                          ParametersImt.mcl)

    def connect_ue_to_bs(self):
        """
        Link the UE randomly to a BS to which the path coupling loss is within
        the minimum coupling loss plus the HO margin.
        """
        for ue in range(self.ue.num_stations):
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
        """
        Calculates interference that IMT system generates on other system
        """

        self.coupling_loss_ue_sat = np.array(np.transpose(
                                self.calculate_coupling_loss(self.ue, self.system,
                                            self.propagation_system)).tolist()[0])
        self.coupling_loss_bs_sat = np.array(np.transpose(
                                self.calculate_coupling_loss(self.bs, self.system,
                                            self.propagation_system)).tolist()[0])

        ue_bandwidth = self.num_rb_per_ue * self.param.rb_bandwidth
        #bs_bandwidth = self.num_rb_per_bs * self.param.rb_bandwidth

        # applying a bandwidth scaling factor since UE transmits on a portion
        # of the satellite's bandwidth
        interference_ue = self.ue.tx_power - self.coupling_loss_ue_sat \
                            + 10*math.log10(ue_bandwidth/self.param_system.sat_bandwidth)

        # assume BS transmits with full power (no power control) in the whole bandwidth
        interference_bs = self.param.bs_tx_power - self.coupling_loss_bs_sat

        self.system.rx_interference = 10*math.log10( \
                        math.pow(10, 0.1*self.system.rx_interference)
                        + np.sum(np.power(10, 0.1*interference_ue)) \
                        + np.sum(np.power(10, 0.1*interference_bs)))

        self.system.thermal_noise = \
            10*math.log10(self.param_system.BOLTZMANN_CONSTANT* \
                          self.param_system.sat_noise_temperature) + \
                          10*np.log10(self.param_system.sat_bandwidth * 1e6)

        self.system.total_interference = \
            10*np.log10(np.power(10, 0.1*self.system.rx_interference) + \
                        np.power(10, 0.1*self.system.thermal_noise))

        self.system.inr = self.system.rx_interference / self.system.thermal_noise


    def collect_results(self):
        self.results.add_coupling_loss_dl( \
            np.reshape(self.coupling_loss, self.ue.num_stations*self.bs.num_stations).tolist())
        self.results.add_coupling_loss_bs_sat(self.coupling_loss_bs_sat.tolist())
        self.results.add_coupling_loss_ue_sat(self.coupling_loss_ue_sat.tolist())

        self.results.add_inr([self.system.inr.tolist()])

        # NOT COLLECTING POWER STATISTICS FOR THE MOMENT
#        for bs in range(self.bs.num_stations):
#            self.results.add_tx_power_dl(self.ue.tx_power[bs])
#
#        # select the active stations
#        ids = np.where(self.bs.active)
#
#        self.results.add_sinr_dl(self.bs.sinr[ids].tolist())
#        self.results.add_snr_dl(self.bs.snr[ids].tolist())

