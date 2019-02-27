# -*- coding: utf-8 -*-
"""
Created on Fri Apr  7 17:02:35 2017

@author: edgar
"""

import numpy as np
import math

from sharc.simulation import Simulation
from sharc.parameters.parameters import Parameters
from sharc.station_factory import StationFactory
from sharc.support.enumerations import StationType

from sharc.propagation.propagation_factory import PropagationFactory

class SimulationUplink(Simulation):
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
        #self.plot_scenario()

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
            #self.add_external_interference()
            #self.recalculate_sinr()
            #self.calculate_imt_degradation()
            pass
        else:
            # Execute this piece of code if IMT generates interference into
            # the other system
            self.calculate_sinr()
            self.calculate_external_interference()
            #self.calculate_external_degradation()
            pass

        self.collect_results(write_to_file, snapshot_number)


    def power_control(self):
        """
        Apply uplink power control algorithm
        """
        if self.parameters.imt.ue_tx_power_control == "OFF":
            ue_active = np.where(self.ue.active)[0]
            self.ue.tx_power[ue_active] = self.parameters.imt.ue_p_cmax * np.ones(len(ue_active))
        else:
            bs_active = np.where(self.bs.active)[0]
            for bs in bs_active:
                ue = self.link[bs]
                p_cmax = self.parameters.imt.ue_p_cmax
                m_pusch = self.num_rb_per_ue
                p_o_pusch = self.parameters.imt.ue_p_o_pusch
                alpha = self.parameters.imt.ue_alpha
                cl = self.coupling_loss_imt[bs,ue] + self.parameters.imt.bs_ohmic_loss \
                            + self.parameters.imt.ue_ohmic_loss + self.parameters.imt.ue_body_loss
                self.ue.tx_power[ue] = np.minimum(p_cmax, 10*np.log10(m_pusch) + p_o_pusch + alpha*cl)
        if self.adjacent_channel: 
            self.ue_power_diff = self.parameters.imt.ue_p_cmax - self.ue.tx_power


    def calculate_sinr(self):
        """
        Calculates the uplink SINR for each BS.
        """
        # calculate uplink received power for each active BS
        bs_active = np.where(self.bs.active)[0]
        for bs in bs_active:
            ue = self.link[bs]

            self.bs.rx_power[bs] = self.ue.tx_power[ue]  \
                                        - self.parameters.imt.ue_ohmic_loss \
                                        - self.parameters.imt.ue_body_loss \
                                        - self.coupling_loss_imt[bs,ue] - self.parameters.imt.bs_ohmic_loss
            # create a list of BSs that serve the interfering UEs
            bs_interf = [b for b in bs_active if b not in [bs]]

            # calculate intra system interference
            for bi in bs_interf:
                ui = self.link[bi]
                interference = self.ue.tx_power[ui] - self.parameters.imt.ue_ohmic_loss  \
                                - self.parameters.imt.ue_body_loss \
                                - self.coupling_loss_imt[bs,ui] - self.parameters.imt.bs_ohmic_loss
                self.bs.rx_interference[bs] = 10*np.log10( \
                    np.power(10, 0.1*self.bs.rx_interference[bs])
                    + np.power(10, 0.1*interference))

            # calculate N
            self.bs.thermal_noise[bs] = \
                10*np.log10(self.parameters.imt.BOLTZMANN_CONSTANT*self.parameters.imt.noise_temperature*1e3) + \
                10*np.log10(self.bs.bandwidth[bs] * 1e6) + \
                self.bs.noise_figure[bs]

            # calculate I+N
            self.bs.total_interference[bs] = \
                10*np.log10(np.power(10, 0.1*self.bs.rx_interference[bs]) + \
                            np.power(10, 0.1*self.bs.thermal_noise[bs]))

            # calculate SNR and SINR
            self.bs.sinr[bs] = self.bs.rx_power[bs] - self.bs.total_interference[bs]
            self.bs.snr[bs] = self.bs.rx_power[bs] - self.bs.thermal_noise[bs]


    def calculate_sinr_ext(self):
        """
        Calculates the downlink SINR for each UE taking into account the
        interference that is generated by the other system into IMT system.
        """
        self.coupling_loss_imt_system = self.calculate_coupling_loss(self.system,
                                                                     self.bs,
                                                                     self.propagation_system)

        bs_active = np.where(self.bs.active)[0]
        tx_power = self.param_system.tx_power_density + 10*np.log10(self.bs.bandwidth*1e6) + 30
        for bs in bs_active:
            active_beams = [i for i in range(bs*self.parameters.imt.ue_k, (bs+1)*self.parameters.imt.ue_k)]
            self.bs.ext_interference[bs] = tx_power[bs] - self.coupling_loss_imt_system[active_beams] \
                                            - self.parameters.imt.bs_ohmic_loss

            self.bs.sinr_ext[bs] = self.bs.rx_power[bs] \
                - (10*np.log10(np.power(10, 0.1*self.bs.total_interference[bs]) + np.power(10, 0.1*self.bs.ext_interference[bs])))
            self.bs.inr[bs] = self.bs.ext_interference[bs] - self.bs.thermal_noise[bs]


    def calculate_external_interference(self):
        """
        Calculates interference that IMT system generates on other system
        """
        if self.co_channel:
            self.coupling_loss_imt_system = self.calculate_coupling_loss(self.system,
                                                                         self.ue,
                                                                         self.propagation_system)

        if self.adjacent_channel:
              self.coupling_loss_imt_system_adjacent = self.calculate_coupling_loss(self.system,
                                                                       self.ue,
                                                                       self.propagation_system,
                                                                       c_channel=False)

        # applying a bandwidth scaling factor since UE transmits on a portion
        # of the satellite's bandwidth
        # calculate interference only from active UE's
        rx_interference = 0

        bs_active = np.where(self.bs.active)[0]
        for bs in bs_active:
            ue = self.link[bs]

            if self.co_channel:
                if self.overlapping_bandwidth:
                    acs = 0
                else:
                    acs = self.param_system.adjacent_ch_selectivity

                interference_ue = self.ue.tx_power[ue] - self.parameters.imt.ue_ohmic_loss \
                                  - self.parameters.imt.ue_body_loss \
                                  - self.coupling_loss_imt_system[ue]
                weights = self.calculate_bw_weights(self.parameters.imt.bandwidth,
                                                    self.param_system.bandwidth,
                                                    self.parameters.imt.ue_k)
                rx_interference += np.sum(weights*np.power(10, 0.1*interference_ue)) / 10**(acs/10.)

            if self.adjacent_channel:
                oob_power = self.ue.spectral_mask.power_calc(self.param_system.frequency,self.system.bandwidth)\
                            - self.ue_power_diff[ue]
                oob_interference_array = oob_power - self.coupling_loss_imt_system_adjacent[ue] \
                                            + 10*np.log10((self.param_system.bandwidth - self.overlapping_bandwidth)/
                                              self.param_system.bandwidth) \
                                            - self.parameters.imt.ue_body_loss
                rx_interference += np.sum(np.power(10,0.1*oob_interference_array))

        self.system.rx_interference = 10*np.log10(rx_interference)
        # calculate N
        self.system.thermal_noise = \
            10*np.log10(self.param_system.BOLTZMANN_CONSTANT* \
                          self.system.noise_temperature*1e3) + \
                          10*math.log10(self.param_system.bandwidth * 1e6)

        # calculate INR at the system
        self.system.inr = np.array([self.system.rx_interference - self.system.thermal_noise])

        # Calculate PFD at the system
        if self.system.station_type is StationType.RAS:
            self.system.pfd = 10*np.log10(10**(self.system.rx_interference/10)/self.system.antenna[0].effective_area)


    def collect_results(self, write_to_file: bool, snapshot_number: int):
        if not self.parameters.imt.interfered_with and np.any(self.bs.active):
            self.results.system_inr.extend(self.system.inr.tolist())
            self.results.system_ul_interf_power.extend([self.system.rx_interference])
            if self.system.station_type is StationType.RAS:
                self.results.system_pfd.extend([self.system.pfd])

        bs_active = np.where(self.bs.active)[0]
        for bs in bs_active:
            ue = self.link[bs]
            self.results.imt_path_loss.extend(self.path_loss_imt[bs,ue])
            self.results.imt_coupling_loss.extend(self.coupling_loss_imt[bs,ue])

            self.results.imt_bs_antenna_gain.extend(self.imt_bs_antenna_gain[bs,ue])
            self.results.imt_ue_antenna_gain.extend(self.imt_ue_antenna_gain[bs,ue])

            tput = self.calculate_imt_tput(self.bs.sinr[bs],
                                           self.parameters.imt.ul_sinr_min,
                                           self.parameters.imt.ul_sinr_max,
                                           self.parameters.imt.ul_attenuation_factor)
            self.results.imt_ul_tput.extend(tput.tolist())

            if self.parameters.imt.interfered_with:
                tput_ext = self.calculate_imt_tput(self.bs.sinr_ext[bs],
                                                      self.parameters.imt.ul_sinr_min,
                                                      self.parameters.imt.ul_sinr_max,
                                                      self.parameters.imt.ul_attenuation_factor)
                self.results.imt_ul_tput_ext.extend(tput_ext.tolist())
                self.results.imt_ul_sinr_ext.extend(self.bs.sinr_ext[bs].tolist())
                self.results.imt_ul_inr.extend(self.bs.inr[bs].tolist())

                active_beams = [i for i in range(bs*self.parameters.imt.ue_k, (bs+1)*self.parameters.imt.ue_k)]
                self.results.system_imt_antenna_gain.extend(self.system_imt_antenna_gain[0,active_beams])
                self.results.imt_system_antenna_gain.extend(self.imt_system_antenna_gain[0,active_beams])
                self.results.imt_system_path_loss.extend(self.imt_system_path_loss[0,active_beams])
                if self.param_system.channel_model == "HDFSS":
                    self.results.imt_system_build_entry_loss.extend(self.imt_system_build_entry_loss[:,bs])
                    self.results.imt_system_diffraction_loss.extend(self.imt_system_diffraction_loss[:,bs])
            else:
                self.results.system_imt_antenna_gain.extend(self.system_imt_antenna_gain[0,ue])
                self.results.imt_system_antenna_gain.extend(self.imt_system_antenna_gain[0,ue])
                self.results.imt_system_path_loss.extend(self.imt_system_path_loss[0,ue])
                if self.param_system.channel_model == "HDFSS":
                    self.results.imt_system_build_entry_loss.extend(self.imt_system_build_entry_loss[:,ue])
                    self.results.imt_system_diffraction_loss.extend(self.imt_system_diffraction_loss[:,ue])

            self.results.imt_ul_tx_power.extend(self.ue.tx_power[ue].tolist())
            imt_ul_tx_power_density = 10*np.log10(np.power(10, 0.1*self.ue.tx_power[ue])/(self.num_rb_per_ue*self.parameters.imt.rb_bandwidth*1e6))
            self.results.imt_ul_tx_power_density.extend(imt_ul_tx_power_density.tolist())
            self.results.imt_ul_sinr.extend(self.bs.sinr[bs].tolist())
            self.results.imt_ul_snr.extend(self.bs.snr[bs].tolist())

        if write_to_file:
            self.results.write_files(snapshot_number)
            self.notify_observers(source=__name__, results=self.results)



