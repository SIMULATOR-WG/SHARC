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

class SimulationDownlink(Simulation):
    """
    Implements the flowchart of simulation downlink method
    """

    def __init__(self, parameters: Parameters):
        super().__init__(parameters)

        
    def snapshot(self, *args, **kwargs):
        write_to_file = kwargs["write_to_file"]
        snapshot_number = kwargs["snapshot_number"]
        
        # In case of hotspots, base stations coordinates have to be calculated
        # on every snapshot. Anyway, let topology decide whether to calculate
        # or not
        self.topology.calculate_coordinates()
        
        # Create the base stations (remember that it takes into account the
        # network load factor)
        self.bs = StationFactory.generate_imt_base_stations(self.parameters.imt,
                                                            self.parameters.antenna_imt,
                                                            self.topology)      
        
        # Create the other system (FSS, HAPS, etc...)
        self.system = StationFactory.generate_system(self.parameters)

        # Create IMT user equipments
        self.ue = StationFactory.generate_imt_ue(self.parameters.imt,
                                                 self.parameters.antenna_imt,
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
        
        if self.parameters.imt.interfered_with:
            # Execute this piece of code if the other system generates 
            # interference into IMT
            self.calculate_sinr()
            self.calculate_sinr_ext()
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


        
    def finalize(self, *args, **kwargs):
        self.notify_observers(source=__name__, results=self.results)


    def power_control(self):
        """
        Apply downlink power control algorithm
        """
        # Currently, the maximum transmit power of the base station is equaly
        # divided among the selected UEs
        tx_power = self.parameters.imt.bs_conducted_power + self.bs_power_gain \
                    - self.parameters.imt.bs_feed_loss - 10*math.log10(self.parameters.imt.ue_k) 
        # calculate tansmit powers to have a structure such as
        # {bs_1: [pwr_1, pwr_2,...], ...}, where bs_1 is the base station id,
        # pwr_1 is the transmit power from bs_1 to ue_1, pwr_2 is the transmit
        # power from bs_1 to ue_2, etc
        bs_active = np.where(self.bs.active)[0]
        self.bs.tx_power = dict([(bs, tx_power*np.ones(self.parameters.imt.ue_k)) for bs in bs_active])

        
    def calculate_sinr(self):
        """
        Calculates the downlink SINR for each UE. 
        """        
        bs_active = np.where(self.bs.active)[0]
        for bs in bs_active:
            ue = self.link[bs]
            self.ue.rx_power[ue] = self.bs.tx_power[bs] - self.coupling_loss_imt[bs,ue] \
                                     - self.parameters.imt.ue_body_loss - self.parameters.imt.ue_feed_loss \
                                     - self.parameters.imt.bs_feed_loss

            # create a list with base stations that generate interference in ue_list
            bs_interf = [b for b in bs_active if b not in [bs]]

            # calculate intra system interference
            for bi in bs_interf:
                interference = self.bs.tx_power[bi] - self.coupling_loss_imt[bi,ue] \
                                 - self.parameters.imt.ue_body_loss - self.parameters.imt.ue_feed_loss \
                                 - self.parameters.imt.bs_feed_loss
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
        self.coupling_loss_imt_system = self.calculate_coupling_loss(self.system, 
                                                                     self.ue,
                                                                     self.propagation_system)       
        
        # applying a bandwidth scaling factor since UE transmits on a portion
        # of the satellite's bandwidth
        # calculate interference only to active UE's
        ue = np.where(self.ue.active)[0]
        tx_power = self.param_system.tx_power_density + 10*np.log10(self.ue.bandwidth[ue]*1e6) + 30
        self.ue.ext_interference[ue] = tx_power - self.coupling_loss_imt_system[ue] \
                            - self.parameters.imt.ue_body_loss - self.parameters.imt.ue_feed_loss

        self.ue.sinr_ext[ue] = self.ue.rx_power[ue] \
            - (10*np.log10(np.power(10, 0.1*self.ue.total_interference[ue]) + np.power(10, 0.1*self.ue.ext_interference[ue])))
        self.ue.inr[ue] = self.ue.ext_interference[ue] - self.ue.thermal_noise[ue]
        
        
    def calculate_external_interference(self):
        """
        Calculates interference that IMT system generates on other system
        """

        self.coupling_loss_imt_system = self.calculate_coupling_loss(self.system, 
                                                                     self.bs,
                                                                     self.propagation_system)

        # applying a bandwidth scaling factor since UE transmits on a portion
        # of the satellite's bandwidth
        # calculate interference only from active UE's
        bs_active = np.where(self.bs.active)[0]
        for bs in bs_active:
            active_beams = [i for i in range(bs*self.parameters.imt.ue_k, (bs+1)*self.parameters.imt.ue_k)]
            interference = self.bs.tx_power[bs] - self.coupling_loss_imt_system[active_beams] \
                                + 10*np.log10(self.bs.bandwidth[bs]/self.param_system.bandwidth)
            weights = self.calculate_bw_weights(self.parameters.imt.bandwidth, 
                                                self.param_system.bandwidth,
                                                self.parameters.imt.ue_k)
            self.system.rx_interference = 10*math.log10( \
                    math.pow(10, 0.1*self.system.rx_interference) + np.sum(weights*np.power(10, 0.1*interference)))

        # calculate N
        self.system.thermal_noise = \
            10*math.log10(self.param_system.BOLTZMANN_CONSTANT* \
                          self.param_system.noise_temperature*1e3) + \
                          10*math.log10(self.param_system.bandwidth * 1e6)

        # calculate INR at the system
        self.system.inr = np.array([self.system.rx_interference - self.system.thermal_noise])
        
        
    def collect_results(self, write_to_file: bool, snapshot_number: int):
        if not self.parameters.imt.interfered_with:
            self.results.system_inr.extend(self.system.inr.tolist())
            self.results.system_inr_scaled.extend([self.system.inr + 10*math.log10(self.param_system.inr_scaling)])
        
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
                
                self.results.system_imt_antenna_gain.extend(self.system_imt_antenna_gain[0,ue])
                self.results.imt_system_antenna_gain.extend(self.imt_system_antenna_gain[0,ue])                
            else:
                active_beams = [i for i in range(bs*self.parameters.imt.ue_k, (bs+1)*self.parameters.imt.ue_k)]
                self.results.system_imt_antenna_gain.extend(self.system_imt_antenna_gain[0,active_beams])
                self.results.imt_system_antenna_gain.extend(self.imt_system_antenna_gain[0,active_beams])

            self.results.imt_dl_tx_power.extend(self.bs.tx_power[bs].tolist())
            #imt_dl_tx_power_density = 10*np.log10(np.power(10, 0.1*self.bs.tx_power[bs])/(self.num_rb_per_ue*self.parameters.imt.rb_bandwidth*1e6))
            #self.results.imt_dl_tx_power_density.extend(imt_dl_tx_power_density.tolist())
            self.results.imt_dl_sinr.extend(self.ue.sinr[ue].tolist())
            self.results.imt_dl_snr.extend(self.ue.snr[ue].tolist())
            
        if write_to_file:
            self.results.write_files(snapshot_number)
            self.notify_observers(source=__name__, results=self.results)

