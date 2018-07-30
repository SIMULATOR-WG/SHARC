# -*- coding: utf-8 -*-
"""
Created on Wed Jan 11 19:04:03 2017

@author: edgar
"""

from abc import ABC, abstractmethod
from sharc.support.observable import Observable

import numpy as np
import math
import sys
import matplotlib.pyplot as plt

from sharc.support.enumerations import StationType
from sharc.topology.topology_factory import TopologyFactory
from sharc.parameters.parameters import Parameters
from sharc.propagation.propagation import Propagation
from sharc.station_manager import StationManager
from sharc.results import Results
from sharc.propagation.propagation_factory import PropagationFactory


class Simulation(ABC, Observable):

    def __init__(self, parameters: Parameters, parameter_file: str):
        ABC.__init__(self)
        Observable.__init__(self)

        self.parameters = parameters
        self.parameters_filename = parameter_file

        if self.parameters.general.system == "FSS_SS":
            self.param_system = self.parameters.fss_ss
        elif self.parameters.general.system == "FSS_ES":
            self.param_system = self.parameters.fss_es
        elif self.parameters.general.system == "FS":
            self.param_system = self.parameters.fs
        elif self.parameters.general.system == "HAPS":
            self.param_system = self.parameters.haps
        elif self.parameters.general.system == "RNS":
            self.param_system = self.parameters.rns
        elif self.parameters.general.system == "RAS":
            self.param_system = self.parameters.ras

        self.co_channel = self.parameters.general.enable_cochannel
        self.adjacent_channel = self.parameters.general.enable_adjacent_channel

        self.topology = TopologyFactory.createTopology(self.parameters)

        self.bs_power_gain = 0
        self.ue_power_gain = 0

        self.imt_bs_antenna_gain = list()
        self.imt_ue_antenna_gain = list()
        self.system_imt_antenna_gain = list()
        self.imt_system_antenna_gain = list()

        self.path_loss_imt = np.empty(0)
        self.coupling_loss_imt = np.empty(0)
        self.coupling_loss_imt_system = np.empty(0)
        self.coupling_loss_imt_system_adjacent = np.empty(0)

        self.bs_to_ue_phi = np.empty(0)
        self.bs_to_ue_theta = np.empty(0)
        self.bs_to_ue_beam_rbs = np.empty(0)

        self.ue = np.empty(0)
        self.bs = np.empty(0)
        self.system = np.empty(0)

        self.link = dict()

        self.num_rb_per_bs = 0
        self.num_rb_per_ue = 0

        self.results = None

        imt_min_freq = self.parameters.imt.frequency - self.parameters.imt.bandwidth / 2
        imt_max_freq = self.parameters.imt.frequency + self.parameters.imt.bandwidth / 2
        system_min_freq = self.param_system.frequency - self.param_system.bandwidth / 2
        system_max_freq = self.param_system.frequency + self.param_system.bandwidth / 2

        max_min_freq = np.maximum(imt_min_freq, system_min_freq)
        min_max_freq = np.minimum(imt_max_freq, system_max_freq)

        self.overlapping_bandwidth = min_max_freq - max_min_freq
        if self.overlapping_bandwidth < 0:
            self.overlapping_bandwidth = 0

        if (self.overlapping_bandwidth == self.param_system.bandwidth and
            not self.parameters.imt.interfered_with) or \
           (self.overlapping_bandwidth == self.parameters.imt.bandwidth and
            self.parameters.imt.interfered_with):

            self.adjacent_channel = False

        self.propagation_imt = None
        self.propagation_system = None

    def add_observer_list(self, observers: list):
        for o in observers:
            self.add_observer(o)

    def initialize(self, *args, **kwargs):
        """
        This method is executed only once to initialize the simulation variables.
        """

        self.topology.calculate_coordinates()
        num_bs = self.topology.num_base_stations
        num_ue = num_bs*self.parameters.imt.ue_k*self.parameters.imt.ue_k_m

        self.bs_power_gain = 10*math.log10(self.parameters.antenna_imt.bs_tx_n_rows*
                                           self.parameters.antenna_imt.bs_tx_n_columns)
        self.ue_power_gain = 10*math.log10(self.parameters.antenna_imt.ue_tx_n_rows*
                                           self.parameters.antenna_imt.ue_tx_n_columns)
        self.imt_bs_antenna_gain = list()
        self.imt_ue_antenna_gain = list()
        self.path_loss_imt = np.empty([num_bs, num_ue])
        self.coupling_loss_imt = np.empty([num_bs, num_ue])
        self.coupling_loss_imt_system = np.empty(num_ue)

        self.bs_to_ue_phi = np.empty([num_bs, num_ue])
        self.bs_to_ue_theta = np.empty([num_bs, num_ue])
        self.bs_to_ue_beam_rbs = -1.0*np.ones(num_ue, dtype=int)

        self.ue = np.empty(num_ue)
        self.bs = np.empty(num_bs)
        self.system = np.empty(1)

        # this attribute indicates the list of UE's that are connected to each
        # base station. The position the the list indicates the resource block
        # group that is allocated to the given UE
        self.link = dict([(bs,list()) for bs in range(num_bs)])

        # calculates the number of RB per BS
        self.num_rb_per_bs = math.trunc((1-self.parameters.imt.guard_band_ratio)* \
                            self.parameters.imt.bandwidth /self.parameters.imt.rb_bandwidth)
        # calculates the number of RB per UE on a given BS
        self.num_rb_per_ue = math.trunc(self.num_rb_per_bs/self.parameters.imt.ue_k)

        self.results = Results(self.parameters_filename, self.parameters.general.overwrite_output)

    def finalize(self, *args, **kwargs):
        """
        Finalizes the simulation (collect final results, etc...)
        """
        snapshot_number = kwargs["snapshot_number"]
        self.results.write_files(snapshot_number)

    def calculate_coupling_loss(self,
                                station_a: StationManager,
                                station_b: StationManager,
                                propagation: Propagation,
                                c_channel = True) -> np.array:
        """
        Calculates the path coupling loss from each station_a to all station_b.
        Result is returned as a numpy array with dimensions num_a x num_b
        TODO: calculate coupling loss between activa stations only
        """
        # Calculate distance from transmitters to receivers. The result is a
        # num_bs x num_ue array
        d_2D = station_a.get_distance_to(station_b)
        d_3D = station_a.get_3d_distance_to(station_b)

        if self.parameters.imt.interfered_with:
            freq = self.param_system.frequency
        else:
            freq = self.parameters.imt.frequency

        if station_a.station_type is StationType.FSS_SS or \
           station_a.station_type is StationType.HAPS or \
           station_a.station_type is StationType.RNS:
            elevation_angles = station_b.get_elevation_angle(station_a, self.param_system)
        elif station_a.station_type is StationType.IMT_BS and \
             station_b.station_type is StationType.IMT_UE and \
             self.parameters.imt.topology == "INDOOR":
            elevation_angles = np.transpose(station_b.get_elevation(station_a))
        elif station_a.station_type is StationType.FSS_ES:
            elevation_angles = station_b.get_elevation(station_a)
        else:
            elevation_angles = None

        if station_a.station_type is StationType.FSS_SS or \
           station_a.station_type is StationType.FSS_ES or \
           station_a.station_type is StationType.HAPS or \
           station_a.station_type is StationType.FS or \
           station_a.station_type is StationType.RNS or \
           station_a.station_type is StationType.RAS:

            if station_b.station_type is StationType.IMT_UE:
                # define antenna gains
                gain_a = self.calculate_gains(station_a, station_b)
                gain_b = np.transpose(self.calculate_gains(station_b, station_a, c_channel))
                sectors_in_node=1

            else:
                # define antenna gains
                gain_a = np.repeat(self.calculate_gains(station_a, station_b), self.parameters.imt.ue_k, 1)
                gain_b = np.transpose(self.calculate_gains(station_b, station_a, c_channel))
                sectors_in_node = self.parameters.imt.ue_k

            if self.parameters.imt.interfered_with:
                earth_to_space = False
                single_entry = True
            else:
                earth_to_space = True
                single_entry = False

            if station_a.station_type is StationType.FSS_SS or \
               station_a.station_type is StationType.HAPS or \
               station_a.station_type is StationType.RNS:
                path_loss = propagation.get_loss(distance_3D=d_3D,
                                             frequency=freq*np.ones(d_3D.shape),
                                             indoor_stations=np.tile(station_b.indoor, (station_a.num_stations, 1)),
                                             elevation=elevation_angles, sat_params = self.param_system,
                                             earth_to_space = earth_to_space, earth_station_antenna_gain=gain_b,
                                             single_entry=single_entry, number_of_sectors=sectors_in_node)
            else:
                path_loss = propagation.get_loss(distance_3D=d_3D,
                                             frequency=freq*np.ones(d_3D.shape),
                                             indoor_stations=np.tile(station_b.indoor, (station_a.num_stations, 1)),
                                             elevation=elevation_angles, es_params=self.param_system,
                                             tx_gain = gain_a, rx_gain = gain_b, number_of_sectors=sectors_in_node)

            self.system_imt_antenna_gain = gain_a
            self.imt_system_antenna_gain = gain_b
        else:
            path_loss = propagation.get_loss(distance_3D=d_3D,
                                             distance_2D=d_2D,
                                             frequency=self.parameters.imt.frequency*np.ones(d_2D.shape),
                                             indoor_stations=np.tile(station_b.indoor, (station_a.num_stations, 1)),
                                             bs_height=station_a.height,
                                             ue_height=station_b.height,
                                             elevation=elevation_angles,
                                             shadowing=self.parameters.imt.shadowing,
                                             line_of_sight_prob=self.parameters.imt.line_of_sight_prob)
            # define antenna gains
            gain_a = self.calculate_gains(station_a, station_b)
            gain_b = np.transpose(self.calculate_gains(station_b, station_a))

            # collect IMT BS and UE antenna gain samples
            self.path_loss_imt = path_loss
            self.imt_bs_antenna_gain = gain_a
            self.imt_ue_antenna_gain = gain_b

        # calculate coupling loss
        coupling_loss = np.squeeze(path_loss - gain_a - gain_b)

        return coupling_loss

    def connect_ue_to_bs(self):
        """
        Link the UE's to the serving BS. It is assumed that each group of K*M
        user equipments are distributed and pointed to a certain base station
        according to the decisions taken at TG 5/1 meeting
        """
        num_ue_per_bs = self.parameters.imt.ue_k*self.parameters.imt.ue_k_m
        bs_active = np.where(self.bs.active)[0]
        for bs in bs_active:
            ue_list = [i for i in range(bs*num_ue_per_bs, bs*num_ue_per_bs + num_ue_per_bs)]
            self.link[bs] = ue_list

    def select_ue(self, random_number_gen: np.random.RandomState):
        """
        Select K UEs randomly from all the UEs linked to one BS as “chosen”
        UEs. These K “chosen” UEs will be scheduled during this snapshot.
        """
        self.bs_to_ue_phi, self.bs_to_ue_theta = \
            self.bs.get_pointing_vector_to(self.ue)

        bs_active = np.where(self.bs.active)[0]
        for bs in bs_active:
            # select K UE's among the ones that are connected to BS
            random_number_gen.shuffle(self.link[bs])
            K = self.parameters.imt.ue_k
            del self.link[bs][K:]
            # Activate the selected UE's and create beams
            if self.bs.active[bs]:
                self.ue.active[self.link[bs]] = np.ones(K, dtype=bool)
                for ue in self.link[bs]:
                    # add beam to BS antennas
                    self.bs.antenna[bs].add_beam(self.bs_to_ue_phi[bs,ue],
                                             self.bs_to_ue_theta[bs,ue])
                    # add beam to UE antennas
                    self.ue.antenna[ue].add_beam(self.bs_to_ue_phi[bs,ue] - 180,
                                             180 - self.bs_to_ue_theta[bs,ue])
                    # set beam resource block group
                    self.bs_to_ue_beam_rbs[ue] = len(self.bs.antenna[bs].beams_list) - 1


    def scheduler(self):
        """
        This scheduler divides the available resource blocks among UE's for
        a given BS
        """
        bs_active = np.where(self.bs.active)[0]
        for bs in bs_active:
            ue = self.link[bs]
            self.bs.bandwidth[bs] = self.num_rb_per_ue*self.parameters.imt.rb_bandwidth
            self.ue.bandwidth[ue] = self.num_rb_per_ue*self.parameters.imt.rb_bandwidth

    def calculate_gains(self,
                        station_1: StationManager,
                        station_2: StationManager,
                        c_channel = True) -> np.array:
        """
        Calculates the gains of antennas in station_1 in the direction of
        station_2
        """
        phi, theta = station_1.get_pointing_vector_to(station_2)

        station_1_active = np.where(station_1.active)[0]
        station_2_active = np.where(station_2.active)[0]

        # Initialize variables (phi, theta, beams_idx)
        if(station_1.station_type is StationType.IMT_BS):
            if(station_2.station_type is StationType.IMT_UE):
                beams_idx = self.bs_to_ue_beam_rbs[station_2_active]
            elif(station_2.station_type is StationType.FSS_SS or \
                 station_2.station_type is StationType.FSS_ES or \
                 station_2.station_type is StationType.HAPS or \
                 station_2.station_type is StationType.FS or \
                 station_2.station_type is StationType.RNS or \
                 station_2.station_type is StationType.RAS):
                phi = np.repeat(phi,self.parameters.imt.ue_k,0)
                theta = np.repeat(theta,self.parameters.imt.ue_k,0)
                beams_idx = np.tile(np.arange(self.parameters.imt.ue_k),self.bs.num_stations)

        elif(station_1.station_type is StationType.IMT_UE):
            beams_idx = np.zeros(len(station_2_active),dtype=int)

        elif(station_1.station_type is StationType.FSS_SS or \
             station_1.station_type is StationType.FSS_ES or \
             station_1.station_type is StationType.HAPS or \
             station_1.station_type is StationType.FS or \
             station_1.station_type is StationType.RNS or \
             station_1.station_type is StationType.RAS):
            beams_idx = np.zeros(len(station_2_active),dtype=int)

        # Calculate gains
        gains = np.zeros(phi.shape)
        if (station_1.station_type is StationType.IMT_BS and station_2.station_type is StationType.FSS_SS) or \
           (station_1.station_type is StationType.IMT_BS and station_2.station_type is StationType.FSS_ES) or \
           (station_1.station_type is StationType.IMT_BS and station_2.station_type is StationType.HAPS) or \
           (station_1.station_type is StationType.IMT_BS and station_2.station_type is StationType.FS) or \
           (station_1.station_type is StationType.IMT_BS and station_2.station_type is StationType.RNS) or \
           (station_1.station_type is StationType.IMT_BS and station_2.station_type is StationType.RAS):
            for k in station_1_active:
                for b in range(k*self.parameters.imt.ue_k,(k+1)*self.parameters.imt.ue_k):
                    gains[b,station_2_active] = station_1.antenna[k].calculate_gain(phi_vec=phi[b,station_2_active],
                                                                            theta_vec=theta[b,station_2_active],
                                                                            beams_l=np.array([beams_idx[b]]),
                                                                            co_channel=c_channel)

        elif (station_1.station_type is StationType.IMT_UE and station_2.station_type is StationType.FSS_SS) or \
             (station_1.station_type is StationType.IMT_UE and station_2.station_type is StationType.FSS_ES) or \
             (station_1.station_type is StationType.IMT_UE and station_2.station_type is StationType.HAPS) or \
             (station_1.station_type is StationType.IMT_UE and station_2.station_type is StationType.FS) or \
             (station_1.station_type is StationType.IMT_UE and station_2.station_type is StationType.RNS) or \
             (station_1.station_type is StationType.IMT_UE and station_2.station_type is StationType.RAS):
               for k in station_1_active:
                   gains[k,station_2_active] = station_1.antenna[k].calculate_gain(phi_vec=phi[k,station_2_active],
                                                                            theta_vec=theta[k,station_2_active],
                                                                            beams_l=beams_idx,
                                                                            co_channel=c_channel)

        elif station_1.station_type is StationType.RNS:
            gains[0,station_2_active] = station_1.antenna[0].calculate_gain(phi_vec = phi[0,station_2_active],
                                                                            theta_vec = theta[0,station_2_active])

        elif station_1.station_type is StationType.FSS_SS or \
             station_1.station_type is StationType.FSS_ES or \
             station_1.station_type is StationType.HAPS or \
             station_1.station_type is StationType.FS or \
             station_1.station_type is StationType.RAS:

            off_axis_angle = station_1.get_off_axis_angle(station_2)
            distance = station_1.get_distance_to(station_2)
            theta = np.degrees(np.arctan((station_1.height - station_2.height)/distance)) + station_1.elevation
            gains[0,station_2_active] = station_1.antenna[0].calculate_gain(off_axis_angle_vec=off_axis_angle[0,station_2_active],
                                                                            theta_vec=theta[0,station_2_active])
        else: # for IMT <-> IMT
            for k in station_1_active:
                gains[k,station_2_active] = station_1.antenna[k].calculate_gain(phi_vec=phi[k,station_2_active],
                                                                            theta_vec=theta[k,station_2_active],
                                                                            beams_l=beams_idx)

        return gains

    def calculate_imt_tput(self,
                           sinr: np.array,
                           sinr_min: float,
                           sinr_max: float,
                           attenuation_factor: float) -> np.array:
        tput_min = 0
        tput_max = attenuation_factor*math.log2(1+math.pow(10, 0.1*sinr_max))

        tput = attenuation_factor*np.log2(1+np.power(10, 0.1*sinr))

        id_min = np.where(sinr < sinr_min)[0]
        id_max = np.where(sinr > sinr_max)[0]

        if len(id_min) > 0:
            tput[id_min] = tput_min
        if len(id_max) > 0:
            tput[id_max] = tput_max

        return tput

    def calculate_bw_weights(self, bw_imt: float, bw_sys: float, ue_k: int) -> np.array:
        """
        Calculates the weight that each resource block group of IMT base stations
        will have when estimating the interference to other systems based on
        the bandwidths of both systems.

        Parameters
        ----------
            bw_imt : bandwidth of IMT system
            bw_sys : bandwidth of other system
            ue_k : number of UE's allocated to each IMT base station; it also
                corresponds to the number of resource block groups

        Returns
        -------
            K-dimentional array of weights
        """

        if bw_imt <= bw_sys:
            weights = np.ones(ue_k)

        elif bw_imt > bw_sys:
            weights = np.zeros(ue_k)

            bw_per_rbg = bw_imt / ue_k

            # number of resource block groups that will have weight equal to 1
            rb_ones = math.floor( bw_sys / bw_per_rbg )

            # weight of the rbg that will generate partial interference
            rb_partial = np.mod( bw_sys, bw_per_rbg ) / bw_per_rbg

            # assign value to weight array
            weights[:rb_ones] = 1
            weights[rb_ones] = rb_partial

        return weights

    def plot_scenario(self):
        fig = plt.figure(figsize=(8,8), facecolor='w', edgecolor='k')
        ax = fig.gca()

        # Plot network topology
        self.topology.plot(ax)

        # Plot user equipments
        ax.scatter(self.ue.x, self.ue.y, color='r', edgecolor="w", linewidth=0.5, label="UE")

        # Plot UE's azimuth
        d = 0.1 * self.topology.cell_radius
        for i in range(len(self.ue.x)):
            plt.plot([self.ue.x[i], self.ue.x[i] + d*math.cos(math.radians(self.ue.azimuth[i]))],
                     [self.ue.y[i], self.ue.y[i] + d*math.sin(math.radians(self.ue.azimuth[i]))],
                     'r-')

        plt.axis('image')
        plt.title("Simulation scenario")
        plt.xlabel("x-coordinate [m]")
        plt.ylabel("y-coordinate [m]")
        plt.legend(loc="upper left", scatterpoints=1)
        plt.tight_layout()
        plt.show()

        if self.parameters.imt.topology == "INDOOR":
            fig = plt.figure(figsize=(8,8), facecolor='w', edgecolor='k')
            ax = fig.gca()

            # Plot network topology
            self.topology.plot(ax,top_view=False)

            # Plot user equipments
            ax.scatter(self.ue.x, self.ue.height, color='r', edgecolor="w", linewidth=0.5, label="UE")

            plt.title("Simulation scenario: side view")
            plt.xlabel("x-coordinate [m]")
            plt.ylabel("z-coordinate [m]")
            plt.legend(loc="upper left", scatterpoints=1)
            plt.tight_layout()
            plt.show()
        
#        sys.exit(0)

    @abstractmethod
    def snapshot(self, *args, **kwargs):
        """
        Performs a single snapshot.
        """
        pass

    @abstractmethod
    def power_control(self):
        """
        Apply downlink power control algorithm
        """

    @abstractmethod
    def collect_results(self, *args, **kwargs):
        """
        Collects results.
        """
        pass
