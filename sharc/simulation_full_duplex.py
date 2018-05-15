# -*- coding: utf-8 -*-
"""
Created on Tue Sep 12 17:21:06 2017

@author: Calil
"""

import numpy as np
import math

from sharc.simulation import Simulation
from sharc.parameters.parameters import Parameters
from sharc.station_factory import StationFactory

class SimulationFullDuplex(Simulation):
    """
    Implements the full duplex simulation
    """

    def __init__(self, parameters: Parameters):
        super().__init__(parameters)
        self.coupling_loss_imt_bs_bs = np.empty(0)
        self.coupling_loss_imt_ue_ue = np.empty(0)
        self.coupling_loss_imt_bs_system = np.empty(0)
        self.coupling_loss_imt_ue_system = np.empty(0)
        self.system_ul_inr = np.empty(0)
        self.system_dl_inr = np.empty(0)

        
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

        # Create IMT user equipments
        self.ue = StationFactory.generate_imt_ue(self.parameters.imt,
                                                 self.parameters.antenna_imt,
                                                 self.topology)
        
        # Create the other system (FSS, HAPS, etc...)
        self.system = StationFactory.generate_system(self.parameters)
        
        #self.plot_scenario()
        
        self.connect_ue_to_bs()
        self.select_ue()
        
        # Calculate coupling loss after beams are created
        self.coupling_loss_imt = self.calculate_coupling_loss(self.bs, 
                                                              self.ue,
                                                              self.propagation_imt)
        
        # UE to UE coupling loss
        self.coupling_loss_imt_ue_ue = self.calculate_coupling_loss(self.ue,
                                                                    self.ue,
                                                                    self.propagation_imt)
        
        # BS to BS coupling loss
        self.coupling_loss_imt_bs_bs = self.calculate_coupling_loss(self.bs,
                                                                    self.bs,
                                                                    self.propagation_imt)
        
        
        # Scheduler which divides the band equally among BSs and UEs
        self.scheduler()
        
        # Stations power control
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


    def power_control(self):
        """
        Apply downling and uplink power control algorithm
        """
        # Downlink Power Control:
        # Currently, the maximum transmit power of the base station is equaly
        # divided among the selected UEs
        tx_power = self.parameters.imt.bs_conducted_power + self.bs_power_gain \
                    - self.parameters.imt.bs_ohmic_loss - 10*math.log10(self.parameters.imt.ue_k) 
        # calculate tansmit powers to have a structure such as
        # {bs_1: [pwr_1, pwr_2,...], ...}, where bs_1 is the base station id,
        # pwr_1 is the transmit power from bs_1 to ue_1, pwr_2 is the transmit
        # power from bs_1 to ue_2, etc
        bs_active = np.where(self.bs.active)[0]
        self.bs.tx_power = dict([(bs, tx_power*np.ones(self.parameters.imt.ue_k)) for bs in bs_active])
        
        # Uplink power control:
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
                alpha = self.parameters.imt.ue_alfa
                cl = self.coupling_loss_imt[bs,ue] + self.parameters.imt.bs_ohmic_loss \
                            + self.parameters.imt.ue_ohmic_loss + self.parameters.imt.ue_body_loss
                self.ue.tx_power[ue] = np.minimum(p_cmax, 10*np.log10(m_pusch) + p_o_pusch + alpha*cl)
    
        
    def calculate_sinr(self):
        """
        Calculates the downlink and uplink SINR for each UE and BS.
        Self-interference is considered
        """    
        ### Downlink SNIR
        bs_active = np.where(self.bs.active)[0]
        for bs in bs_active:
            ue = self.link[bs]
            self.ue.rx_power[ue] = self.bs.tx_power[bs] - self.coupling_loss_imt[bs,ue] \
                                     - self.parameters.imt.ue_body_loss - self.parameters.imt.ue_ohmic_loss

            # create a list with base stations that generate interference in ue_list
            bs_interf = [b for b in bs_active if b not in [bs]]

            #  Internal interference
            for bi in bs_interf:
                #  Interference from BSs
                interference_bs = self.bs.tx_power[bi] - self.coupling_loss_imt[bi,ue] \
                                 - self.parameters.imt.ue_body_loss - self.parameters.imt.ue_ohmic_loss
                            
                ue_interf = self.link[bi]
                # Interference from UEs
                interference_ue = self.ue.tx_power[ue_interf] - self.coupling_loss_imt_ue_ue[ue_interf,ue] \
                                 - 2*self.parameters.imt.ue_body_loss - self.parameters.imt.ue_ohmic_loss
           
                self.ue.rx_interference[ue] = 10*np.log10( \
                    np.power(10, 0.1*self.ue.rx_interference[ue]) + \
                    np.power(10, 0.1*interference_bs) + \
                    np.power(10, 0.1*interference_ue))
                
            self.ue.self_interference[ue] = self.ue.tx_power[ue] - self.ue.sic[ue]

        self.ue.thermal_noise = \
            10*math.log10(self.parameters.imt.BOLTZMANN_CONSTANT*self.parameters.imt.noise_temperature*1e3) + \
            10*np.log10(self.ue.bandwidth * 1e6) + \
            self.ue.noise_figure
            
        self.ue.total_interference = \
            10*np.log10(np.power(10, 0.1*self.ue.rx_interference) + \
                        np.power(10, 0.1*self.ue.thermal_noise)   + \
                        np.power(10, 0.1*self.ue.self_interference))
            
        self.ue.sinr = self.ue.rx_power - self.ue.total_interference
        self.ue.snr = self.ue.rx_power - self.ue.thermal_noise
        
        ### Uplink SNIR
        bs_active = np.where(self.bs.active)[0]
        for bs in bs_active:
            ue = self.link[bs]
            self.bs.rx_power[bs] = self.ue.tx_power[ue]  \
                                        - self.parameters.imt.ue_ohmic_loss - self.parameters.imt.ue_body_loss \
                                        - self.coupling_loss_imt[bs,ue] - self.parameters.imt.bs_ohmic_loss
            # create a list of BSs that serve the interfering UEs
            bs_interf = [b for b in bs_active if b not in [bs]]

            # calculate intra system interference
            for bi in bs_interf:
                ui = self.link[bi]
                interference_ue = self.ue.tx_power[ui] \
                                - self.parameters.imt.ue_ohmic_loss - self.parameters.imt.ue_body_loss \
                                - self.coupling_loss_imt[bs,ui] - self.parameters.imt.bs_ohmic_loss
                                
                interference_bs = self.bs.tx_power[bi] - self.parameters.imt.bs_ohmic_loss \
                                - self.coupling_loss_imt_bs_bs[bs,np.arange(bi*self.parameters.imt.ue_k,(bi+1)*self.parameters.imt.ue_k)]
                                
                self.bs.rx_interference[bs] = 10*np.log10( \
                    np.power(10, 0.1*self.bs.rx_interference[bs])
                    + np.power(10, 0.1*interference_ue) \
                    + np.power(10, 0.1*interference_bs))
                
            # calculate self interference
            self.bs.self_interference[bs] = self.bs.tx_power[bs] - self.bs.sic[bs]
            
            # calculate N
            self.bs.thermal_noise[bs] = \
                10*np.log10(self.parameters.imt.BOLTZMANN_CONSTANT*self.parameters.imt.noise_temperature*1e3) + \
                10*np.log10(self.bs.bandwidth[bs] * 1e6) + \
                self.bs.noise_figure[bs]
    
            # calculate I+N
            self.bs.total_interference[bs] = \
                10*np.log10(np.power(10, 0.1*self.bs.rx_interference[bs]) + \
                            np.power(10, 0.1*self.bs.thermal_noise[bs])   + \
                            np.power(10, 0.1*self.bs.self_interference[bs]))
                
            # calculate SNR and SINR
            self.bs.sinr[bs] = self.bs.rx_power[bs] - self.bs.total_interference[bs]
            self.bs.snr[bs] = self.bs.rx_power[bs] - self.bs.thermal_noise[bs]

        
    def calculate_sinr_ext(self):
        """
        Calculates the SINR and INR for each UE and BS taking into
        account the interference that is generated by the other system into 
        IMT system.
        """
        ### Downlink
        self.coupling_loss_imt_ue_system = self.calculate_coupling_loss(self.system, 
                                                                        self.ue,
                                                                        self.propagation_system)       
        
        # applying a bandwidth scaling factor since UE transmits on a portion
        # of the satellite's bandwidth
        # calculate interference only to active UE's
        ue = np.where(self.ue.active)[0]
        tx_power = self.param_system.tx_power_density + 10*np.log10(self.ue.bandwidth[ue]*1e6) + 30
        self.ue.ext_interference[ue] = tx_power - self.coupling_loss_imt_ue_system[ue] \
                            - self.parameters.imt.ue_body_loss - self.parameters.imt.ue_ohmic_loss

        self.ue.sinr_ext[ue] = self.ue.rx_power[ue] \
            - (10*np.log10(np.power(10, 0.1*self.ue.total_interference[ue]) + np.power(10, 0.1*self.ue.ext_interference[ue])))
        self.ue.inr[ue] = self.ue.ext_interference[ue] - self.ue.thermal_noise[ue]
        
        ### Uplink
        self.coupling_loss_imt_bs_system = self.calculate_coupling_loss(self.system, 
                                                                     self.bs,
                                                                     self.propagation_system)       
        
        # applying a bandwidth scaling factor since UE transmits on a portion
        # of the satellite's bandwidth
        # calculate interference only to active UE's
        bs_active = np.where(self.bs.active)[0]
        tx_power = self.param_system.tx_power_density + 10*np.log10(self.bs.bandwidth*1e6) + 30
        for bs in bs_active:
            active_beams = [i for i in range(bs*self.parameters.imt.ue_k, (bs+1)*self.parameters.imt.ue_k)]
            self.bs.ext_interference[bs] = tx_power[bs] - self.coupling_loss_imt_bs_system[active_beams] \
                                            - self.parameters.imt.bs_ohmic_loss

            self.bs.sinr_ext[bs] = self.bs.rx_power[bs] \
                - (10*np.log10(np.power(10, 0.1*self.bs.total_interference[bs]) + np.power(10, 0.1*self.bs.ext_interference[bs])))
            self.bs.inr[bs] = self.bs.ext_interference[bs] - self.bs.thermal_noise[bs]
        
        
    def calculate_external_interference(self):
        """
        Calculates interference that IMT system generates on other system
        """
        ### Downlink & Uplink
        self.coupling_loss_imt_bs_system = self.calculate_coupling_loss(self.system, 
                                                                        self.bs,
                                                                        self.propagation_system)
        self.coupling_loss_imt_ue_system = self.calculate_coupling_loss(self.system, 
                                                                        self.ue,
                                                                        self.propagation_system)
        
        # calculate N
        self.system.thermal_noise = \
            10*math.log10(self.param_system.BOLTZMANN_CONSTANT* \
                          self.param_system.noise_temperature*1e3) + \
                          10*math.log10(self.param_system.bandwidth * 1e6)


        # applying a bandwidth scaling factor since UE transmits on a portion
        # of the satellite's bandwidth
        # calculate interference only from active UE's
        bs_active = np.where(self.bs.active)[0]
        for bs in bs_active:
            active_beams = [i for i in range(bs*self.parameters.imt.ue_k, (bs+1)*self.parameters.imt.ue_k)]
            interference_bs = self.bs.tx_power[bs] - self.coupling_loss_imt_bs_system[active_beams] \
                                + 10*np.log10(self.bs.bandwidth[bs]/self.param_system.bandwidth)
                                
            self.system.rx_interference = 10*math.log10( \
                    math.pow(10, 0.1*self.system.rx_interference) + \
                    np.sum(np.power(10, 0.1*interference_bs)))
        
        if not self.parameters.imt.interfered_with:
            self.system_dl_inr = np.array([self.system.rx_interference - self.system.thermal_noise])
        
        # UE interference
        ue_active = np.where(self.ue.active)[0]
        interference_ue = self.ue.tx_power[ue_active] \
                            - self.parameters.imt.ue_ohmic_loss - self.parameters.imt.ue_body_loss \
                            - self.coupling_loss_imt_ue_system[ue_active] \
                            + 10*np.log10(self.ue.bandwidth[ue_active]/self.param_system.bandwidth)
                            
        self.system.rx_interference = 10*np.log10(np.power(10, 0.1*self.system.rx_interference) + \
                                                  np.sum(np.power(10, 0.1*interference_ue)))
        
        if not self.parameters.imt.interfered_with:
            self.system_ul_inr = np.array(interference_ue - self.system.thermal_noise)

        # calculate INR at the system
        self.system.inr = np.array([self.system.rx_interference - self.system.thermal_noise])
        
        
    def collect_results(self, write_to_file: bool, snapshot_number: int):
        if not self.parameters.imt.interfered_with:
            self.results.system_inr.extend(self.system.inr.tolist())
            self.results.system_inr_scaled.extend([self.system.inr + 10*math.log10(self.param_system.inr_scaling)])
            self.results.system_ul_inr_scaled.extend(self.system_ul_inr + 10*math.log10(self.param_system.inr_scaling))
            self.results.system_dl_inr_scaled.extend([self.system_dl_inr + 10*math.log10(self.param_system.inr_scaling)])
        
        bs_active = np.where(self.bs.active)[0]
        for bs in bs_active:
            ue = self.link[bs]
            self.results.imt_path_loss.extend(self.path_loss_imt[bs,ue])
            self.results.imt_coupling_loss.extend(self.coupling_loss_imt[bs,ue])
            
            if not self.parameters.general.suppress_large_results:
                self.results.imt_coupling_loss_all.extend(self.coupling_loss_imt[bs,:])
            
                bs_bs_pl = self.path_loss_imt_bs_bs[bs, ~np.isnan(self.coupling_loss_imt_bs_bs[bs,:])]
                self.results.imt_bs_bs_path_loss.extend(bs_bs_pl)
            
                ue_ue_pl = np.ravel(self.path_loss_imt_ue_ue[ue])
                ue_ue_pl = ue_ue_pl[~np.isnan(ue_ue_pl)]
                self.results.imt_ue_ue_path_loss.extend(ue_ue_pl)
            
                bs_bs_cl = self.coupling_loss_imt_bs_bs[bs, ~np.isnan(self.coupling_loss_imt_bs_bs[bs,:])]
                self.results.imt_coupling_loss_bs_bs.extend(bs_bs_cl)
            
                ue_ue_cl = np.ravel(self.coupling_loss_imt_ue_ue[ue])
                ue_ue_cl = ue_ue_cl[~np.isnan(ue_ue_cl)]
                self.results.imt_coupling_loss_ue_ue.extend(ue_ue_cl)
            
                bs_bs_ag = self.imt_bs_bs_antenna_gain[bs, ~np.isnan(self.imt_bs_bs_antenna_gain[bs,:])]
                self.results.imt_bs_bs_antenna_gain.extend(bs_bs_ag)
            
                ue_ue_ag = np.ravel(self.imt_ue_ue_antenna_gain[ue])
                ue_ue_ag = ue_ue_ag[~np.isnan(ue_ue_ag)]
                self.results.imt_ue_ue_antenna_gain.extend(ue_ue_ag)
            
            self.results.imt_bs_antenna_gain.extend(self.imt_bs_antenna_gain[bs,ue])
            self.results.imt_ue_antenna_gain.extend(self.imt_ue_antenna_gain[bs,ue])
            
            
            tput = self.calculate_imt_tput(self.ue.sinr[ue],
                                           self.parameters.imt.dl_sinr_min,
                                           self.parameters.imt.dl_sinr_max,
                                           self.parameters.imt.dl_attenuation_factor)
            tput = tput*self.ue.bandwidth[ue]
            self.results.imt_dl_tput.extend(tput.tolist())
            
            total_tput = np.sum(tput)
            
            tput = self.calculate_imt_tput(self.bs.sinr[bs],
                                           self.parameters.imt.ul_sinr_min,
                                           self.parameters.imt.ul_sinr_max,
                                           self.parameters.imt.ul_attenuation_factor)
            tput = tput*self.bs.bandwidth[bs]
            self.results.imt_ul_tput.extend(tput.tolist())
            
            total_tput = np.sum(tput) + total_tput
            self.results.imt_total_tput.extend([total_tput])

            if self.parameters.imt.interfered_with:
                tput_ext = self.calculate_imt_tput(self.ue.sinr_ext[ue],
                                                   self.parameters.imt.dl_sinr_min,
                                                   self.parameters.imt.dl_sinr_max,
                                                   self.parameters.imt.dl_attenuation_factor)
                tput_ext = tput*self.ue.bandwidth[ue]
                self.results.imt_dl_tput_ext.extend(tput_ext.tolist()) 
                self.results.imt_dl_sinr_ext.extend(self.ue.sinr_ext[ue].tolist())
                self.results.imt_dl_inr.extend(self.ue.inr[ue].tolist())
                
                tput_ext = self.calculate_imt_tput(self.bs.sinr_ext[bs],
                                                      self.parameters.imt.ul_sinr_min,
                                                      self.parameters.imt.ul_sinr_max,
                                                      self.parameters.imt.ul_attenuation_factor)
                tput_ext = tput*self.bs.bandwidth[bs]
                self.results.imt_ul_tput_ext.extend([tput_ext])  
                self.results.imt_ul_sinr_ext.extend(self.bs.sinr_ext[bs].tolist())
                self.results.imt_ul_inr.extend(self.bs.inr[bs].tolist())
                
            self.results.system_imt_ue_antenna_gain.extend(self.system_imt_ue_antenna_gain[0,ue])
            self.results.imt_ue_system_antenna_gain.extend(self.imt_ue_system_antenna_gain[0,ue])                
                
            active_beams = [i for i in range(bs*self.parameters.imt.ue_k, (bs+1)*self.parameters.imt.ue_k)]
            self.results.system_imt_bs_antenna_gain.extend(self.system_imt_bs_antenna_gain[0,active_beams])
            self.results.imt_bs_system_antenna_gain.extend(self.imt_bs_system_antenna_gain[0,active_beams])
            
            self.results.system_ul_coupling_loss.extend(self.coupling_loss_imt_ue_system[ue])
            self.results.system_dl_coupling_loss.extend([self.coupling_loss_imt_bs_system[bs]])

            self.results.imt_dl_tx_power.extend(self.bs.tx_power[bs].tolist())
            self.results.imt_dl_rx_power.extend(self.ue.rx_power[ue].tolist())
            self.results.imt_dl_sinr.extend(self.ue.sinr[ue].tolist())
            self.results.imt_dl_snr.extend(self.ue.snr[ue].tolist())
            self.results.imt_dl_ue_interf.extend(self.ue.total_interference[ue].tolist())
            
            self.results.imt_ul_tx_power.extend(self.ue.tx_power[ue].tolist())
            self.results.imt_ul_rx_power.extend(self.bs.rx_power[bs].tolist())
            self.results.imt_ul_sinr.extend(self.bs.sinr[bs].tolist())
            self.results.imt_ul_snr.extend(self.bs.snr[bs].tolist())
            self.results.imt_ul_bs_interf.extend(self.bs.total_interference[bs].tolist())
            
        if write_to_file:
            self.results.write_files(snapshot_number)
            self.notify_observers(source=__name__, results=self.results)

