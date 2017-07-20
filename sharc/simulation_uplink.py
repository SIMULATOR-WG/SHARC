# -*- coding: utf-8 -*-
"""
Created on Fri Apr  7 17:02:35 2017

@author: edgar
"""

import numpy as np
import math

from sharc.simulation import Simulation
from sharc.parameters.parameters_imt import ParametersImt
from sharc.parameters.parameters_antenna_imt import ParametersAntennaImt
from sharc.parameters.parameters_fss import ParametersFss
from sharc.station_factory import StationFactory

class SimulationUplink(Simulation):
    """
    Implements the flowchart of simulation downlink method
    """

    def __init__(self, param_imt: ParametersImt, param_system: ParametersFss, param_ant: ParametersAntennaImt):
        super().__init__(param_imt, param_system, param_ant)

        
    def snapshot(self, *args, **kwargs):
        write_to_file = kwargs["write_to_file"]
        snapshot_number = kwargs["snapshot_number"]
        
        # In case of hotspots, base stations coordinates have to be calculated
        # on every snapshot. Anyway, let topology decide whether to calculate
        # or not
        self.topology.calculate_coordinates()
        
        # Create the base stations (remember that it takes into account the
        # network load factor)
        self.bs = StationFactory.generate_imt_base_stations(self.param_imt,
                                                            self.param_imt_antenna,
                                                            self.topology)      
        
        # Create the other system (FSS, HAPS, etc...)
        # Currently it supports only FSS space station
        self.system = StationFactory.generate_fss_space_stations(self.param_system)

        # Create IMT user equipments
        self.ue = StationFactory.generate_imt_ue(self.param_imt,
                                                 self.param_imt_antenna,
                                                 self.topology)
        #self.plot_scenario()
        
        self.connect_ue_to_bs()
        self.select_ue()
        
        # Calculate coupling loss after beams are created
        self.coupling_loss_imt = self.calculate_coupling_loss(self.bs, 
                                                              self.ue,
                                                              self.propagation_imt)
        self.scheduler()
        self.power_control()
        
        if self.param_imt.interfered_with:
            # Execute this piece of code if the other system generates 
            # interference into IMT
            #self.calculate_sinr()
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
        if self.param_imt.ue_tx_power_control == "OFF":
            ue_active = np.where(self.ue.active)[0]
            self.ue.tx_power[ue_active] = self.param_imt.ue_p_cmax * np.ones(len(ue_active))
        else:
            bs_active = np.where(self.bs.active)[0]
            for bs in bs_active:
                ue = self.link[bs]
                p_cmax = self.param_imt.ue_p_cmax
                m_pusch = self.num_rb_per_ue
                p_o_pusch = self.param_imt.ue_p_o_pusch
                alpha = self.param_imt.ue_alfa
                cl = self.coupling_loss_imt[bs,ue] + self.ue_power_gain
                self.ue.tx_power[ue] = np.minimum(p_cmax, 10*np.log10(m_pusch) + p_o_pusch + alpha*cl)


    def calculate_sinr(self):
        """
        Calculates the uplink SINR for each UE. 
        """
        # calculate uplink received power for each active BS
        bs_active = np.where(self.bs.active)[0]
        for bs in bs_active:
            ue = self.link[bs]
            self.bs.rx_power[bs] = self.ue.tx_power[ue]  \
                                        - self.param_imt.ue_feed_loss - self.param_imt.ue_body_loss \
                                        - self.coupling_loss_imt[bs,ue] - self.param_imt.bs_feed_loss
            # create a list of BSs that serve the interfering UEs
            bs_interf = [b for b in bs_active if b not in [bs]]

            # calculate intra system interference
            for bi in bs_interf:
                ui = self.link[bi]
                interference = self.ue.tx_power[ui] \
                                - self.param_imt.ue_feed_loss - self.param_imt.ue_body_loss \
                                - self.coupling_loss_imt[bs,ui] - self.param_imt.bs_feed_loss
                self.bs.rx_interference[bs] = 10*np.log10( \
                    np.power(10, 0.1*self.bs.rx_interference[bs])
                    + np.power(10, 0.1*interference))

            # calculate N
            self.bs.thermal_noise[bs] = \
                10*np.log10(self.param_imt.BOLTZMANN_CONSTANT*self.param_imt.noise_temperature*1e3) + \
                10*np.log10(self.bs.bandwidth[bs] * 1e6) + \
                self.bs.noise_figure[bs]
    
            # calculate I+N
            self.bs.total_interference[bs] = \
                10*np.log10(np.power(10, 0.1*self.bs.rx_interference[bs]) + \
                            np.power(10, 0.1*self.bs.thermal_noise[bs]))
                
            # calculate SNR and SINR
            self.bs.sinr[bs] = self.bs.rx_power[bs] - self.bs.total_interference[bs]
            self.bs.snr[bs] = self.bs.rx_power[bs] - self.bs.thermal_noise[bs]


    def calculate_external_interference(self):
        """
        Calculates interference that IMT system generates on other system
        """

        self.coupling_loss_imt_system = self.calculate_coupling_loss(self.system, 
                                                                     self.ue,
                                                                     self.propagation_system)

        ue_bandwidth = self.num_rb_per_ue * self.param_imt.rb_bandwidth

        # applying a bandwidth scaling factor since UE transmits on a portion
        # of the satellite's bandwidth
        # calculate interference only from active UE's
        ue_active = np.where(self.ue.active)[0]
        interference_ue = self.ue.tx_power[ue_active] - self.coupling_loss_imt_system[ue_active] \
                            + 10*math.log10(ue_bandwidth/self.param_system.bandwidth) \
                            - self.param_imt.ue_body_loss - self.param_imt.ue_feed_loss

        # calculate the aggregate interference on system
        self.system.rx_interference = 10*math.log10(np.sum(np.power(10, 0.1*interference_ue)))

        # calculate N
        self.system.thermal_noise = \
            10*math.log10(self.param_system.BOLTZMANN_CONSTANT* \
                          self.param_system.sat_noise_temperature) + \
                          10*math.log10(self.param_system.bandwidth * 1e6)

        # calculate INR at the system
        self.system.inr = self.system.rx_interference - self.system.thermal_noise


    def collect_results(self, write_to_file: bool, snapshot_number: int):
        self.results.system_inr.extend([self.system.inr])
        
        bs_active = np.where(self.bs.active)[0]
        for bs in bs_active:
            ue = self.link[bs]
            self.results.imt_path_loss.extend(self.path_loss_imt[bs,ue])
            self.results.imt_coupling_loss.extend(self.coupling_loss_imt[bs,ue])
            self.results.imt_bs_antenna_gain.extend(self.imt_bs_antenna_gain[bs,ue])
            self.results.imt_ue_antenna_gain.extend(self.imt_ue_antenna_gain[bs,ue])
            tput = self.calculate_imt_ul_tput(self.bs.sinr[bs])
            self.results.imt_ul_tput.extend(tput.tolist())
            self.results.imt_ul_tx_power.extend(self.ue.tx_power[ue].tolist())
            imt_ul_tx_power_density = 10*np.log10(np.power(10, 0.1*self.ue.tx_power[ue])/(self.num_rb_per_ue*self.param_imt.rb_bandwidth*1e6))
            self.results.imt_ul_tx_power_density.extend(imt_ul_tx_power_density.tolist())
            self.results.imt_ul_sinr.extend(self.bs.sinr[bs].tolist())
            self.results.imt_ul_snr.extend(self.bs.snr[bs].tolist())
            
        if write_to_file:
            self.results.write_files(snapshot_number)
            self.notify_observers(source=__name__, results=self.results)


    
        
    def calculate_imt_ul_tput(self, sinr: np.array) -> np.array:
        tput_min = 0
        tput_max = self.param_imt.ul_attenuation_factor*math.log2(1+math.pow(10, 0.1*self.param_imt.ul_sinr_max))
        
        tput = self.param_imt.ul_attenuation_factor*np.log2(1+np.power(10, 0.1*sinr))
        
        id_min = np.where(sinr<self.param_imt.ul_sinr_min)[0]
        id_max = np.where(sinr>self.param_imt.ul_sinr_max)[0]

        if len(id_min) > 0:
            tput[id_min] = tput_min
        if len(id_max) > 0:
            tput[id_max] = tput_max

        return tput
        
        

        