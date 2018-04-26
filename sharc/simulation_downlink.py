# -*- coding: utf-8 -*-
"""
Created on Wed Jan 11 19:06:41 2017

@author: edgar
"""

import numpy as np
import math

from sharc.simulation import Simulation
from sharc.parameters.parameters import Parameters
from sharc.station_factory import StationFactory
from sharc.support.enumerations import StationType

from sharc.propagation.propagation_factory import PropagationFactory


class SimulationDownlink(Simulation):
    """
    Implements the flowchart of simulation downlink method
    """

    def __init__(self, parameters: Parameters, parameter_file: str):
        super().__init__(parameters, parameter_file)

    def snapshot(self, *args, **kwargs):
        write_to_file = kwargs["write_to_file"]
        snapshot_number = kwargs["snapshot_number"]
        seed = kwargs["seed"]

        random_number_gen = np.random.RandomState(seed)

        self.propagation_imt = PropagationFactory.create_propagation(self.parameters.imt.channel_model, self.parameters,
                                                                    random_number_gen)
        self.propagation_system = PropagationFactory.create_propagation(self.param_system.channel_model, self.parameters,
                                                                       random_number_gen)

        # In case of hotspots, base stations coordinates have to be calculated
        # on every snapshot. Anyway, let topology decide whether to calculate
        # or not
        self.topology.calculate_coordinates(random_number_gen)

        # Create the base stations (remember that it takes into account the
        # network load factor)
        self.bs = StationFactory.generate_imt_base_stations(self.parameters.imt,
                                                            self.parameters.antenna_imt,
                                                            self.topology, random_number_gen)

        # Create the other system (FSS, HAPS, etc...)
        self.system = StationFactory.generate_system(self.parameters, self.topology, random_number_gen)

        # Create IMT user equipments
        self.ue = StationFactory.generate_imt_ue(self.parameters.imt,
                                                 self.parameters.antenna_imt,
                                                 self.topology, random_number_gen)

        self.connect_ue_to_bs()
        self.select_ue(random_number_gen)

        # Calculate coupling loss after beams are created
        self.coupling_loss_imt = self.calculate_coupling_loss(self.bs,
                                                              self.ue,
                                                              self.propagation_imt)
        self.scheduler()
        self.power_control()

        if self.parameters.imt.interfered_with:
            # Execute this piece of code if the other system generates
            # interference into IMT
            self.calculate_sinr()
            self.calculate_sinr_ext()
            pass
        else:
            # Execute this piece of code if IMT generates interference into
            # the other system
            self.calculate_sinr()
            self.calculate_external_interference()
            pass

        self.collect_results(write_to_file, snapshot_number)

    def finalize(self, *args, **kwargs):
        self.notify_observers(source=__name__, results=self.results)

    def power_control(self):
        """
        Apply downlink power control algorithm
        """
        # Currently, the maximum transmit power of the base station is equaly
        # divided among the selected UEs
        total_power = self.parameters.imt.bs_conducted_power \
                      + self.bs_power_gain
        tx_power = total_power - 10 * math.log10(self.parameters.imt.ue_k)
        # calculate transmit powers to have a structure such as
        # {bs_1: [pwr_1, pwr_2,...], ...}, where bs_1 is the base station id,
        # pwr_1 is the transmit power from bs_1 to ue_1, pwr_2 is the transmit
        # power from bs_1 to ue_2, etc
        bs_active = np.where(self.bs.active)[0]
        self.bs.tx_power = dict([(bs, tx_power*np.ones(self.parameters.imt.ue_k)) for bs in bs_active])

        # Update the spectral mask
        if self.adjacent_channel:
            self.bs.spectral_mask.set_mask(power = total_power)

    def calculate_sinr(self):
        """
        Calculates the downlink SINR for each UE.
        """
        bs_active = np.where(self.bs.active)[0]
        for bs in bs_active:
            ue = self.link[bs]
            self.ue.rx_power[ue] = self.bs.tx_power[bs] - self.parameters.imt.bs_ohmic_loss \
                                       - self.coupling_loss_imt[bs,ue] \
                                       - self.parameters.imt.ue_body_loss \
                                       - self.parameters.imt.ue_ohmic_loss

            # create a list with base stations that generate interference in ue_list
            bs_interf = [b for b in bs_active if b not in [bs]]

            # calculate intra system interference
            for bi in bs_interf:
                interference = self.bs.tx_power[bi] - self.parameters.imt.bs_ohmic_loss \
                                 - self.coupling_loss_imt[bi,ue] \
                                 - self.parameters.imt.ue_body_loss - self.parameters.imt.ue_ohmic_loss

                self.ue.rx_interference[ue] = 10*np.log10( \
                    np.power(10, 0.1*self.ue.rx_interference[ue]) + np.power(10, 0.1*interference))

        self.ue.thermal_noise = \
            10*math.log10(self.parameters.imt.BOLTZMANN_CONSTANT*self.parameters.imt.noise_temperature*1e3) + \
            10*np.log10(self.ue.bandwidth * 1e6) + \
            self.ue.noise_figure

        self.ue.total_interference = \
            10*np.log10(np.power(10, 0.1*self.ue.rx_interference) + \
                        np.power(10, 0.1*self.ue.thermal_noise))

        self.ue.sinr = self.ue.rx_power - self.ue.total_interference
        self.ue.snr = self.ue.rx_power - self.ue.thermal_noise

    def calculate_sinr_ext(self):
        """
        Calculates the downlink SINR and INR for each UE taking into account the
        interference that is generated by the other system into IMT system.
        """
        self.coupling_loss_imt_system = \
            self.calculate_coupling_loss_haps_platform(self.system,
                                                       self.ue,
                                                       self.propagation_system,
                                                       c_channel = self.co_channel)

        # applying a bandwidth scaling factor since UE transmits on a portion
        # of the satellite's bandwidth
        # calculate interference only to active UE's
        ue = np.where(self.ue.active)[0]

        tx_power_sys = self.param_system.tx_power_density + 10*np.log10(self.ue.bandwidth[ue[0]]*1e6)
        rx_interference = tx_power_sys - self.coupling_loss_imt_system[:,ue] \
                            - self.parameters.imt.ue_body_loss - self.parameters.imt.ue_ohmic_loss
        self.ue.ext_interference[ue] = 10*np.log10(np.sum(10**(0.1*rx_interference), 0))
        self.ue.sinr_ext[ue] = self.ue.rx_power[ue] \
            - (10*np.log10(np.power(10, 0.1*self.ue.total_interference[ue]) + np.power(10, 0.1*self.ue.ext_interference[ue])))
        self.ue.inr[ue] = self.ue.ext_interference[ue] - self.ue.thermal_noise[ue]
        protection_criteria = -6
        lambda_imt = 299792458/(self.parameters.imt.frequency*1e6)
        self.ue.pfd[ue] = protection_criteria + 10*np.log10(4*np.pi/(lambda_imt**2)) \
                            - self.imt_system_antenna_gain[0,ue] - 174 + self.parameters.imt.ue_noise_figure
        pfd_level = self.param_system.tx_power_density + self.system_imt_antenna_gain[:,ue] \
                    - self.path_loss_imt_system[:,ue] - self.parameters.imt.ue_body_loss - 3
        self.ue.pfd_level[ue] = 10*np.log10(np.sum(10**(0.1*pfd_level), 0))
        self.ue.pfd_interfered[ue] = self.ue.pfd[ue] < self.ue.pfd_level[ue]

    def calculate_external_interference(self):
        """
        Calculates interference that IMT system generates on other system
        """
        polarization_loss = 3

        if self.co_channel:
            self.coupling_loss_imt_system = self.calculate_coupling_loss(self.system,
                                                                     self.bs,
                                                                     self.propagation_system) + polarization_loss

        if self.adjacent_channel:
            self.coupling_loss_imt_system_adjacent = self.calculate_coupling_loss(self.system,
                                                                     self.bs,
                                                                     self.propagation_system,
                                                                     c_channel=False) + polarization_loss

        # applying a bandwidth scaling factor since UE transmits on a portion
        # of the interfered systems bandwidth
        # calculate interference only from active UE's
        rx_interference = 0

        bs_active = np.where(self.bs.active)[0]
        for bs in bs_active:

            active_beams = [i for i in range(bs*self.parameters.imt.ue_k, (bs+1)*self.parameters.imt.ue_k)]

            if self.co_channel:
                if self.overlapping_bandwidth:
                    acs = 0
                else:
                    acs = self.param_system.adjacent_ch_selectivity

                interference = self.bs.tx_power[bs] - self.parameters.imt.bs_ohmic_loss \
                             - self.coupling_loss_imt_system[active_beams]
                weights = self.calculate_bw_weights(self.parameters.imt.bandwidth,
                                                    self.param_system.bandwidth,
                                                    self.parameters.imt.ue_k)

                rx_interference += np.sum(weights*np.power(10, 0.1*interference)) / 10**(acs/10.)

            if self.adjacent_channel:

                oob_power = self.bs.spectral_mask.power_calc(self.param_system.frequency,self.system.bandwidth)

                oob_interference = oob_power - self.coupling_loss_imt_system_adjacent[active_beams[0]] \
                                   + 10*np.log10((self.param_system.bandwidth - self.overlapping_bandwidth)/
                                                 self.param_system.bandwidth)
                                   
                rx_interference += math.pow(10, 0.1*oob_interference)

        self.system.rx_interference = 10*np.log10(rx_interference)
        # calculate N
        self.system.thermal_noise = \
            10*math.log10(self.param_system.BOLTZMANN_CONSTANT* \
                          self.system.noise_temperature*1e3) + \
                          10*math.log10(self.param_system.bandwidth * 1e6)

        # calculate INR at the system
        self.system.inr = np.array([self.system.rx_interference - self.system.thermal_noise])

        # Calculate PFD at the system
        if self.system.station_type is StationType.RAS:
            self.system.pfd = 10*np.log10(10**(self.system.rx_interference/10)/self.system.antenna[0].effective_area)

    def collect_results(self, write_to_file: bool, snapshot_number: int):
        if not self.parameters.imt.interfered_with:
            self.results.system_inr.extend(self.system.inr.tolist())
            self.results.system_inr_scaled.extend([self.system.inr + 10*math.log10(self.param_system.inr_scaling)])
            if self.system.station_type is StationType.RAS:
                self.results.system_pfd.extend([self.system.pfd])
                self.results.system_dl_interf_power.extend([self.system.rx_interference])
        else:
            self.results.imt_dl_pfd_interfered += np.count_nonzero(self.ue.pfd_interfered)
            self.results.imt_dl_pfd_total += np.count_nonzero(self.ue.active)
            

        bs_active = np.where(self.bs.active)[0]
        for bs in bs_active:
            ue = self.link[bs]
            self.results.imt_path_loss.extend(self.path_loss_imt[bs,ue])
            self.results.imt_coupling_loss.extend(self.coupling_loss_imt[bs,ue])

            self.results.imt_bs_antenna_gain.extend(self.imt_bs_antenna_gain[bs,ue])
            self.results.imt_ue_antenna_gain.extend(self.imt_ue_antenna_gain[bs,ue])

            tput = self.calculate_imt_tput(self.ue.sinr[ue],
                                           self.parameters.imt.dl_sinr_min,
                                           self.parameters.imt.dl_sinr_max,
                                           self.parameters.imt.dl_attenuation_factor)
            self.results.imt_dl_tput.extend(tput.tolist())

            if self.parameters.imt.interfered_with:
                tput_ext = self.calculate_imt_tput(self.ue.sinr_ext[ue],
                                                   self.parameters.imt.dl_sinr_min,
                                                   self.parameters.imt.dl_sinr_max,
                                                   self.parameters.imt.dl_attenuation_factor)
                self.results.imt_dl_tput_ext.extend(tput_ext.tolist())
                self.results.imt_dl_sinr_ext.extend(self.ue.sinr_ext[ue].tolist())
                self.results.imt_dl_inr.extend(self.ue.inr[ue].tolist())
                self.results.imt_dl_pfd.extend(self.ue.pfd[ue].tolist())
                self.results.imt_dl_pfd_level.extend(self.ue.pfd_level[ue].tolist())
                self.results.system_imt_antenna_gain.extend(self.system_imt_antenna_gain[:,ue].flatten().tolist())
                self.results.imt_system_antenna_gain.extend(self.imt_system_antenna_gain[:,ue].flatten().tolist())
            else:
                active_beams = [i for i in range(bs*self.parameters.imt.ue_k, (bs+1)*self.parameters.imt.ue_k)]
                self.results.system_imt_antenna_gain.extend(self.system_imt_antenna_gain[0,active_beams])
                self.results.imt_system_antenna_gain.extend(self.imt_system_antenna_gain[0,active_beams])

            self.results.imt_dl_tx_power.extend(self.bs.tx_power[bs].tolist())

            self.results.imt_dl_sinr.extend(self.ue.sinr[ue].tolist())
            self.results.imt_dl_snr.extend(self.ue.snr[ue].tolist())

        if write_to_file:
            self.results.write_files(snapshot_number)
            self.notify_observers(source=__name__, results=self.results)

