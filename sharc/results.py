# -*- coding: utf-8 -*-
"""
Created on Fri May 11 08:17:55 2018

@author: Calil
"""

import abc

class Results(object):
    
    __metaclass__ = abc.ABCMeta
    
    def __init__(self, out_dir):
        self.imt_ul_tx_power_density = list()
        self.imt_ul_tx_power = list()
        self.imt_ul_sinr_ext = list()
        self.imt_ul_sinr = list()
        self.imt_ul_snr = list()
        self.imt_ul_inr = list()
        self.imt_ul_tput_ext = list()
        self.imt_ul_tput = list()
        
        self.imt_total_tput = list()

        self.imt_path_loss = list()
        self.imt_coupling_loss = list()
        self.imt_coupling_loss_all = list()
        self.imt_bs_antenna_gain = list()
        self.imt_ue_antenna_gain = list()
        
        self.imt_coupling_loss_bs_bs = list()
        self.imt_coupling_loss_ue_ue = list()
        
        self.system_imt_antenna_gain = list()
        self.imt_system_antenna_gain = list()
        
        self.system_imt_ue_antenna_gain = list()
        self.system_imt_bs_antenna_gain = list()
        self.imt_ue_system_antenna_gain = list()
        self.imt_bs_system_antenna_gain = list()

        self.imt_dl_tx_power_density = list()
        self.imt_dl_tx_power = list()
        self.imt_dl_sinr_ext = list()
        self.imt_dl_sinr = list()
        self.imt_dl_snr = list()
        self.imt_dl_inr = list()
        self.imt_dl_tput_ext = list()
        self.imt_dl_tput = list()
        
        self.imt_bs_bs_path_loss = list()
        self.imt_ue_ue_path_loss = list()
        
        self.imt_bs_bs_antenna_gain = list()
        self.imt_ue_ue_antenna_gain = list()
        
        self.system_ul_coupling_loss = list()
        self.system_ul_interf_power = list()
        self.system_ul_inr_scaled = list()

        self.system_dl_coupling_loss = list()
        self.system_dl_interf_power = list()
        self.system_dl_inr_scaled = list()
        
        self.imt_dl_rx_power = list()
        self.imt_ul_rx_power = list()
        self.imt_dl_ue_interf = list()
        self.imt_ul_bs_interf = list()

        self.system_inr = list()
        self.system_inr_scaled = list()
        self.output_directory = out_dir
        
        self.plot_list = list()

    @abc.abstractmethod
    def generate_plot_list(self, n_bins):
        pass            
            
    @abc.abstractmethod
    def write_files(self, snapshot_number: int):
        pass

      