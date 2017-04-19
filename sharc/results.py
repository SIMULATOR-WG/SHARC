# -*- coding: utf-8 -*-
"""
Created on Thu Mar 23 08:47:46 2017

@author: edgar
"""

import numpy as np
import os

class Results(object):
    
    def __init__(self):
        self.tx_power_ul = list()
        self.sinr_ul = list()
        self.snr_ul = list()
        self.throughoput_ul = list()
        self.inr = list()
        self.coupling_loss_ul = list()
        self.coupling_loss_ue_sat = list()
        self.coupling_loss_bs_sat = list()
        
    def add_coupling_loss_ul(self, sample):
        self.coupling_loss_ul.extend(sample)
        
    def add_tx_power_ul(self, sample):
        self.tx_power_ul.extend(sample)
        
    def add_sinr_ul(self, sample):
        self.sinr_ul.extend(sample)
        
    def add_snr_ul(self, sample):
        self.snr_ul.extend(sample)        
        
    def add_inr(self, sample):
        self.inr.extend(sample)
        
    def add_coupling_loss_ue_sat(self, sample):
        self.coupling_loss_ue_sat.extend(sample)
        
    def add_coupling_loss_bs_sat(self, sample):
        self.coupling_loss_bs_sat.extend(sample)
        
    def write_files(self):
        #np.savetxt(os.path.join("results", "coupling_loss_ul.txt"), np.array(self.coupling_loss_ul))
        
        np.savetxt(os.path.join("output", "coupling_loss_ue_sat.txt"), 
                   np.array(self.coupling_loss_ue_sat), fmt="%10.5f")
        np.savetxt(os.path.join("output", "coupling_loss_bs_sat.txt"), 
                   np.array(self.coupling_loss_bs_sat), fmt="%10.5f")
        np.savetxt(os.path.join("output", "inr.txt"), 
                   np.array(self.inr), fmt="%10.5f")
        