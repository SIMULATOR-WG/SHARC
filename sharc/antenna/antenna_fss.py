# -*- coding: utf-8 -*-
"""
Created on Thu Apr 13 17:18:59 2017

@author: edgar
"""

from sharc.antenna.antenna import Antenna
from sharc.parameters.parameters_fss import ParametersFss

import numpy as np

class AntennaFss(Antenna):
    """
    Implements the satellite antenna pattern in the fixed-satellite service
    according to Recommendation ITU-R S.672-4
    """
    
    def __init__(self, param: ParametersFss):
        super(Antenna, self).__init__()
        self.peak_gain = param.sat_rx_antenna_gain
        self.l_s = param.sat_rx_antenna_l_s
        
        if self.l_s == -20:
            self.a = 2.58
        elif self.l_s == -25:
            self.a = 2.88
        elif self.l_s == -30:
            self.a = 3.16
        else:
            pass
            # wrong parameter, write a validation test here

        self.b = 6.32
        
        self.psi_0 = param.sat_rx_antenna_3_dB/2
        self.psi_1 = self.psi_0 * np.power(10, (self.peak_gain + self.l_s + 20)/25)
            
        
    def calculate_gain(self, psi: np.array):
        
        p = np.absolute(psi)
        
        num_samples = len(p)
        gain = np.empty(num_samples)
        
        idx_0 = np.where(p < self.psi_0)
        gain[idx_0] = self.peak_gain
        
        idx_1 = np.intersect1d(np.where(self.psi_0 <= p), np.where(p <= self.a * self.psi_0))
        gain[idx_1] = self.peak_gain - 3 * np.power(p/self.psi_0, 2)
            
        idx_2 = np.intersect1d(np.where(self.a * self.psi_0 < p), np.where(p <= self.b * self.psi_0))
        gain[idx_2] = self.peak_gain + self.l_s
        
        idx_3 = np.intersect1d(np.where(self.b * self.psi_0 < p), np.where(p <= self.psi_1))
        gain[idx_3] = self.peak_gain + self.l_s + 20 - 25 * np.log10(p/self.psi_0)

        idx_4 = np.where(self.psi_1 < p)
        gain[idx_4] = 0
        
        self.gain = gain