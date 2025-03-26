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
from sharc.station_manager import StationManager
from sharc.results import Results
from sharc.propagation.propagation_factory import PropagationFactory


class Simulation(ABC, Observable):

    def __init__(self, parameters: Parameters, parameter_file: str):
        ABC.__init__(self)
        Observable.__init__(self)

        self.parameters = parameters
        self.parameters_filename = parameter_file

        if self.parameters.general.system == "METSAT_SS":
            self.param_system = self.parameters.metsat_ss
        elif self.parameters.general.system == "EESS_SS":
            self.param_system = self.parameters.eess_ss
        elif self.parameters.general.system == "SINGLE_EARTH_STATION":
            self.param_system = self.parameters.single_earth_station
        elif self.parameters.general.system == "FSS_SS":
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
        else:
            sys.stderr.write(
                "ERROR\nInvalid system: " +
                self.parameters.general.system,
            )
            sys.exit(1)

        self.wrap_around_enabled = False
        if self.parameters.imt.topology.type == "MACROCELL":
            self.wrap_around_enabled = self.parameters.imt.topology.macrocell.wrap_around \
                                    and self.parameters.imt.topology.macrocell.num_clusters == 1
        if self.parameters.imt.topology.type == "HOTSPOT":
            self.wrap_around_enabled = self.parameters.imt.topology.hotspot.wrap_around \
                                    and self.parameters.imt.topology.hotspot.num_clusters == 1

        self.co_channel = self.parameters.general.enable_cochannel
        self.adjacent_channel = self.parameters.general.enable_adjacent_channel

        self.topology = TopologyFactory.createTopology(self.parameters)

        self.bs_power_gain = 0
        self.ue_power_gain = 0

        self.imt_bs_antenna_gain = list()
        self.imt_ue_antenna_gain = list()
        self.system_imt_antenna_gain = list()
        self.imt_system_antenna_gain = list()
        self.imt_system_path_loss = list()
        self.imt_system_build_entry_loss = list()
        self.imt_system_diffraction_loss = list()

        self.path_loss_imt = np.empty(0)
        self.coupling_loss_imt = np.empty(0)
        self.coupling_loss_imt_system = np.empty(0)
        self.coupling_loss_imt_system_adjacent = np.empty(0)

        self.bs_to_ue_d_2D = np.empty(0)
        self.bs_to_ue_d_3D = np.empty(0)
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

        if (self.overlapping_bandwidth == self.param_system.bandwidth and not self.parameters.imt.interfered_with) or \
                (self.overlapping_bandwidth == self.parameters.imt.bandwidth and self.parameters.imt.interfered_with):

            self.adjacent_channel = False

        if not self.co_channel and not self.adjacent_channel:
            raise ValueError("Both co_channel and adjacent_channel can't be false")

        random_number_gen = np.random.RandomState(self.parameters.general.seed)
        self.propagation_imt = PropagationFactory.create_propagation(
            self.parameters.imt.channel_model,
            self.parameters,
            self.parameters.imt,
            random_number_gen,
        )
        self.propagation_system = PropagationFactory.create_propagation(
            self.param_system.channel_model,
            self.parameters,
            self.param_system,
            random_number_gen,
        )

    def add_observer_list(self, observers: list):
        for o in observers:
            self.add_observer(o)

    def initialize(self, *args, **kwargs):
        """
        This method is executed only once to initialize the simulation variables.
        """

        self.topology.calculate_coordinates()
        num_bs = self.topology.num_base_stations
        num_ue = num_bs * self.parameters.imt.ue.k * self.parameters.imt.ue.k_m

        self.bs_power_gain = 10 * math.log10(
            self.parameters.imt.bs.antenna.n_rows *
            self.parameters.imt.bs.antenna.n_columns,
        )
        self.ue_power_gain = 10 * math.log10(
            self.parameters.imt.ue.antenna.n_rows *
            self.parameters.imt.ue.antenna.n_columns,
        )
        self.imt_bs_antenna_gain = list()
        self.imt_ue_antenna_gain = list()
        self.path_loss_imt = np.empty([num_bs, num_ue])
        self.coupling_loss_imt = np.empty([num_bs, num_ue])
        self.coupling_loss_imt_system = np.empty(num_ue)

        self.bs_to_ue_phi = np.empty([num_bs, num_ue])
        self.bs_to_ue_theta = np.empty([num_bs, num_ue])
        self.bs_to_ue_beam_rbs = -1.0 * np.ones(num_ue, dtype=int)

        self.ue = np.empty(num_ue)
        self.bs = np.empty(num_bs)
        self.system = np.empty(1)

        # this attribute indicates the list of UE's that are connected to each
        # base station. The position the the list indicates the resource block
        # group that is allocated to the given UE
        self.link = dict([(bs, list()) for bs in range(num_bs)])

        # calculates the number of RB per BS
        self.num_rb_per_bs = math.trunc(
            (1 - self.parameters.imt.guard_band_ratio) *
            self.parameters.imt.bandwidth / self.parameters.imt.rb_bandwidth,
        )
        # calculates the number of RB per UE on a given BS
        self.num_rb_per_ue = math.trunc(
            self.num_rb_per_bs / self.parameters.imt.ue.k,
        )

        self.results = Results().prepare_to_write(
            self.parameters_filename,
            self.parameters.general.overwrite_output,
            self.parameters.general.output_dir,
            self.parameters.general.output_dir_prefix,
        )

        if hasattr(self.param_system, "polarization_loss"):
            self.polarization_loss = self.param_system.polarization_loss
        else:
            self.polarization_loss = 3.0

    def finalize(self, *args, **kwargs):
        """
        Finalizes the simulation (collect final results, etc...)
        """
        snapshot_number = kwargs["snapshot_number"]
        self.results.write_files(snapshot_number)

    def calculate_coupling_loss_system_imt(
        self,
        system_station: StationManager,
        imt_station: StationManager,
        is_co_channel=True,
    ) -> np.array:
        """
        Calculates the coupling loss (path loss + antenna gains + other losses) between
        a system station and an IMT station.

        Returns an numpy array with system_station.size X imt_station.size with coupling loss
        values.

        Parameters
        ----------
        system_station : StationManager
            A StationManager object with system stations
        imt_station : StationManager
            A StationManager object with IMT stations
        is_co_channel : bool, optional
            Whether the interference analysis is co-channel or not, by default True

        Returns
        -------
        np.array
            Returns an numpy array with system_station.size X imt_station.size with coupling loss
            values.
        """
        # Set the frequency and other parameters for the propagation model
        if self.parameters.imt.interfered_with:
            freq = self.param_system.frequency
        else:
            freq = self.parameters.imt.frequency

        # Calculate the antenna gains of the IMT station with respect to the system's station
        if imt_station.station_type is StationType.IMT_UE:
            # define antenna gains
            gain_sys_to_imt = self.calculate_gains(system_station, imt_station)
            gain_imt_to_sys = np.transpose(
                self.calculate_gains(
                    imt_station, system_station, is_co_channel,
                ),
            )
            additional_loss = self.parameters.imt.ue.ohmic_loss \
                + self.parameters.imt.ue.body_loss \
                + self.polarization_loss
        elif imt_station.station_type is StationType.IMT_BS:
            # define antenna gains
            # repeat for each BS beam
            gain_sys_to_imt = np.repeat(
                self.calculate_gains(system_station, imt_station),
                self.parameters.imt.ue.k, 1,
            )
            gain_imt_to_sys = np.transpose(
                self.calculate_gains(
                    imt_station, system_station, is_co_channel,
                ),
            )
            additional_loss = self.parameters.imt.bs.ohmic_loss \
                + self.polarization_loss
        else:
            # should never reach this line
            return ValueError(f"Invalid IMT StationType! {imt_station.station_type}")

        # Calculate the path loss based on the propagation model
        path_loss = self.propagation_system.get_loss(
            self.parameters,
            freq,
            system_station,
            imt_station,
            gain_sys_to_imt,
            gain_imt_to_sys,
        )
        # Store antenna gains and path loss samples
        if self.param_system.channel_model == "HDFSS":
            self.imt_system_build_entry_loss = path_loss[1]
            self.imt_system_diffraction_loss = path_loss[2]
            path_loss = path_loss[0]

        if imt_station.station_type is StationType.IMT_UE:
            self.imt_system_path_loss = path_loss
        else:
            # Repeat for each BS beam
            self.imt_system_path_loss = np.repeat(
                path_loss, self.parameters.imt.ue.k, 1,
            )

        self.system_imt_antenna_gain = gain_sys_to_imt
        self.imt_system_antenna_gain = gain_imt_to_sys

        # calculate coupling loss
        coupling_loss = np.squeeze(
            self.imt_system_path_loss - self.system_imt_antenna_gain -
            self.imt_system_antenna_gain,
        ) + additional_loss

        return coupling_loss

    def calculate_intra_imt_coupling_loss(
        self,
        imt_ue_station: StationManager,
        imt_bs_station: StationManager,
    ) -> np.array:
        """
        Calculates the coupling loss (path loss + antenna gains + other losses) between
        a IMT stations (UE and BS).

        Returns an numpy array with imt_bs_station.size X imt_ue_station.size with coupling loss
        values.

        Parameters
        ----------
        system_station : StationManager
            A StationManager object representins IMT_UE stations
        imt_station : StationManager
            A StationManager object representins IMT_BS stations
        is_co_channel : bool, optional
            Whether the interference analysis is co-channel or not, by default True.
            This parameter is ignored. It's keeped to maintein method interface.

        Returns
        -------
        np.array
            Returns an numpy array with imt_bs_station.size X imt_ue_station.size with coupling loss
            values.
        """
        # Calculate the antenna gains

        ant_gain_bs_to_ue = self.calculate_gains(
            imt_bs_station, imt_ue_station,
        )
        ant_gain_ue_to_bs = self.calculate_gains(
            imt_ue_station, imt_bs_station,
        )

        # Calculate the path loss between IMT stations. Primarly used for UL power control.

        # Note on the array dimentions for coupling loss calculations:
        # The function get_loss returns an array station_a x station_b
        path_loss = self.propagation_imt.get_loss(
            self.parameters,
            self.parameters.imt.frequency,
            imt_ue_station,
            imt_bs_station,
            ant_gain_ue_to_bs,
            ant_gain_bs_to_ue,
        )

        # Collect IMT BS and UE antenna gain samples
        self.path_loss_imt = np.transpose(path_loss)
        self.imt_bs_antenna_gain = ant_gain_bs_to_ue
        self.imt_ue_antenna_gain = np.transpose(ant_gain_ue_to_bs)
        additional_loss = self.parameters.imt.bs.ohmic_loss \
            + self.parameters.imt.ue.ohmic_loss \
            + self.parameters.imt.ue.body_loss

        # calculate coupling loss
        coupling_loss = np.squeeze(
            self.path_loss_imt - self.imt_bs_antenna_gain - self.imt_ue_antenna_gain,
        ) + additional_loss

        return coupling_loss

    def connect_ue_to_bs(self):
        """
        Link the UE's to the serving BS. It is assumed that each group of K*M
        user equipments are distributed and pointed to a certain base station
        according to the decisions taken at TG 5/1 meeting
        """
        num_ue_per_bs = self.parameters.imt.ue.k * self.parameters.imt.ue.k_m
        bs_active = np.where(self.bs.active)[0]
        for bs in bs_active:
            ue_list = [
                i for i in range(
                    bs * num_ue_per_bs, bs * num_ue_per_bs + num_ue_per_bs,
                )
            ]
            self.link[bs] = ue_list

    def select_ue(self, random_number_gen: np.random.RandomState):
        """
        Select K UEs randomly from all the UEs linked to one BS as “chosen”
        UEs. These K “chosen” UEs will be scheduled during this snapshot.
        """
        if self.wrap_around_enabled:
            self.bs_to_ue_d_2D, self.bs_to_ue_d_3D, self.bs_to_ue_phi, self.bs_to_ue_theta = \
                self.bs.get_dist_angles_wrap_around(self.ue)
        else:
            self.bs_to_ue_d_2D = self.bs.get_distance_to(self.ue)
            self.bs_to_ue_d_3D = self.bs.get_3d_distance_to(self.ue)
            self.bs_to_ue_phi, self.bs_to_ue_theta = self.bs.get_pointing_vector_to(
                self.ue,
            )

        bs_active = np.where(self.bs.active)[0]
        for bs in bs_active:
            # select K UE's among the ones that are connected to BS
            random_number_gen.shuffle(self.link[bs])
            K = self.parameters.imt.ue.k
            del self.link[bs][K:]
            # Activate the selected UE's and create beams
            if self.bs.active[bs]:
                self.ue.active[self.link[bs]] = np.ones(K, dtype=bool)
                for ue in self.link[bs]:
                    # add beam to BS antennas
                    self.bs.antenna[bs].add_beam(
                        self.bs_to_ue_phi[bs, ue],
                        self.bs_to_ue_theta[bs, ue],
                    )
                    # add beam to UE antennas
                    self.ue.antenna[ue].add_beam(
                        self.bs_to_ue_phi[bs, ue] - 180,
                        180 - self.bs_to_ue_theta[bs, ue],
                    )
                    # set beam resource block group
                    self.bs_to_ue_beam_rbs[ue] = len(
                        self.bs.antenna[bs].beams_list,
                    ) - 1

    def scheduler(self):
        """
        This scheduler divides the available resource blocks among UE's for
        a given BS
        """
        bs_active = np.where(self.bs.active)[0]
        for bs in bs_active:
            ue = self.link[bs]
            self.bs.bandwidth[bs] = self.num_rb_per_ue * \
                self.parameters.imt.rb_bandwidth
            self.ue.bandwidth[ue] = self.num_rb_per_ue * \
                self.parameters.imt.rb_bandwidth

    def calculate_gains(
        self,
        station_1: StationManager,
        station_2: StationManager,
        c_channel=True,
    ) -> np.array:
        """
        Calculates the gains of antennas in station_1 in the direction of
        station_2
        """
        station_1_active = np.where(station_1.active)[0]
        station_2_active = np.where(station_2.active)[0]

        # Initialize variables (phi, theta, beams_idx)
        if (station_1.station_type is StationType.IMT_BS):
            if (station_2.station_type is StationType.IMT_UE):
                phi = self.bs_to_ue_phi
                theta = self.bs_to_ue_theta
                beams_idx = self.bs_to_ue_beam_rbs[station_2_active]
            elif not station_2.is_imt_station():
                phi, theta = station_1.get_pointing_vector_to(station_2)
                phi = np.repeat(phi, self.parameters.imt.ue.k, 0)
                theta = np.repeat(theta, self.parameters.imt.ue.k, 0)
                beams_idx = np.tile(
                    np.arange(self.parameters.imt.ue.k), self.bs.num_stations,
                )

        elif (station_1.station_type is StationType.IMT_UE):
            phi, theta = station_1.get_pointing_vector_to(station_2)
            beams_idx = np.zeros(len(station_2_active), dtype=int)

        elif not station_1.is_imt_station():
            phi, theta = station_1.get_pointing_vector_to(station_2)
            beams_idx = np.zeros(len(station_2_active), dtype=int)

        # Calculate gains
        gains = np.zeros(phi.shape)
        if station_1.station_type is StationType.IMT_BS and not station_2.is_imt_station():
            for k in station_1_active:
                for b in range(k * self.parameters.imt.ue.k, (k + 1) * self.parameters.imt.ue.k):
                    gains[b, station_2_active] = station_1.antenna[k].calculate_gain(
                        phi_vec=phi[b, station_2_active],
                        theta_vec=theta[
                            b,
                            station_2_active,
                        ],
                        beams_l=np.array(
                            [beams_idx[b]],
                        ),
                        co_channel=c_channel,
                    )

        elif station_1.station_type is StationType.IMT_UE and not station_2.is_imt_station():
            for k in station_1_active:
                gains[k, station_2_active] = station_1.antenna[k].calculate_gain(
                    phi_vec=phi[k, station_2_active],
                    theta_vec=theta[
                        k,
                        station_2_active,
                    ],
                    beams_l=beams_idx,
                    co_channel=c_channel,
                )

        elif station_1.station_type is StationType.RNS:
            gains[0, station_2_active] = station_1.antenna[0].calculate_gain(
                phi_vec=phi[0, station_2_active],
                theta_vec=theta[0, station_2_active],
            )

        elif not station_1.is_imt_station():

            off_axis_angle = station_1.get_off_axis_angle(station_2)
            distance = station_1.get_distance_to(station_2)
            theta = np.degrees(
                np.arctan2(
                    (station_1.height - station_2.height), distance,
                ),
            ) + station_1.elevation
            gains[0, station_2_active] = \
                station_1.antenna[0].calculate_gain(
                    off_axis_angle_vec=off_axis_angle[0, station_2_active],
                    theta_vec=theta[0, station_2_active],
            )
        else:  # for IMT <-> IMT
            for k in station_1_active:
                gains[k, station_2_active] = station_1.antenna[k].calculate_gain(
                    phi_vec=phi[k, station_2_active],
                    theta_vec=theta[
                        k,
                        station_2_active,
                    ],
                    beams_l=beams_idx,
                )
        return gains

    def calculate_imt_tput(
        self,
        sinr: np.array,
        sinr_min: float,
        sinr_max: float,
        attenuation_factor: float,
    ) -> np.array:
        tput_min = 0
        tput_max = attenuation_factor * \
            math.log2(1 + math.pow(10, 0.1 * sinr_max))

        tput = attenuation_factor * np.log2(1 + np.power(10, 0.1 * sinr))

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
            rb_ones = math.floor(bw_sys / bw_per_rbg)

            # weight of the rbg that will generate partial interference
            rb_partial = np.mod(bw_sys, bw_per_rbg) / bw_per_rbg

            # assign value to weight array
            weights[:rb_ones] = 1
            weights[rb_ones] = rb_partial

        return weights

    def plot_scenario(self):
        fig = plt.figure(figsize=(8, 8), facecolor='w', edgecolor='k')
        ax = fig.gca()

        # Plot network topology
        self.topology.plot(ax)

        # Plot user equipments
        ax.scatter(
            self.ue.x, self.ue.y, color='r',
            edgecolor="w", linewidth=0.5, label="UE",
        )

#        wedge = Wedge((0, 0), 300, 0, 360, 290, color='b', alpha=0.2, fill=True)
#        ax.add_artist(wedge)

        # Plot UE's azimuth
        d = 0.1 * self.topology.cell_radius
        for i in range(len(self.ue.x)):
            plt.plot(
                [self.ue.x[i], self.ue.x[i] + d *
                    math.cos(math.radians(self.ue.azimuth[i]))],
                [
                    self.ue.y[i], self.ue.y[i] + d *
                    math.sin(math.radians(self.ue.azimuth[i])),
                ],
                'r-',
            )

        plt.axis('image')
        plt.title("Simulation scenario")
        plt.xlabel("x-coordinate [m]")
        plt.ylabel("y-coordinate [m]")
        plt.legend(loc="upper left", scatterpoints=1)
        plt.tight_layout()
        plt.show()

        if self.parameters.imt.topology == "INDOOR":
            fig = plt.figure(figsize=(8, 8), facecolor='w', edgecolor='k')
            ax = fig.gca()

            # Plot network topology
            self.topology.plot(ax, top_view=False)

            # Plot user equipments
            ax.scatter(
                self.ue.x, self.ue.height, color='r',
                edgecolor="w", linewidth=0.5, label="UE",
            )

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
