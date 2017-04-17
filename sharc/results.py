# -*- coding: utf-8 -*-
"""
Created on Thu Mar 23 08:47:46 2017

@author: edgar
"""

import numpy as np
import os

class Results(object):
    
    def __init__(self):
        self.tx_power_dl = list()
        self.sinr_dl = list()
        self.snr_dl = list()
        self.throughoput_dl = list()
        self.inr = list()
        self.coupling_loss_dl = list()
        self.coupling_loss_ue_sat = list()
        self.coupling_loss_bs_sat = list()
        
    def add_coupling_loss_dl(self, sample):
        self.coupling_loss_dl.extend(sample)
        
    def add_tx_power_dl(self, sample):
        self.tx_power_dl.extend(sample)
        
    def add_sinr_dl(self, sample):
        self.sinr_dl.extend(sample)
        
    def add_snr_dl(self, sample):
        self.snr_dl.extend(sample)        
        
    def add_inr(self, sample):
        self.inr.extend(sample)
        
    def add_coupling_loss_ue_sat(self, sample):
        self.coupling_loss_ue_sat.extend(sample)
        
    def add_coupling_loss_bs_sat(self, sample):
        self.coupling_loss_bs_sat.extend(sample)
        
    def write_files(self):
        #np.savetxt(os.path.join("results", "coupling_loss_dl.txt"), np.array(self.coupling_loss_dl))
        np.savetxt(os.path.join("output", "coupling_loss_ue_sat.txt"), 
                   np.array(self.coupling_loss_ue_sat), fmt="%10.5f")
        np.savetxt(os.path.join("output", "coupling_loss_bs_sat.txt"), 
                   np.array(self.coupling_loss_bs_sat), fmt="%10.5f")
        np.savetxt(os.path.join("output", "inr.txt"), 
                   np.array(self.inr), fmt="%10.5f")
        