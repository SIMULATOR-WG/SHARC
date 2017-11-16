# -*- coding: utf-8 -*-
"""
Created on Thu Nov 16 17:34:11 2017

@author: Calil
"""

from sharc.antenna.antenna import Antenna
from sharc.parameters.parameters_ras import ParametersRas

import numpy as np

class AntenaSA509(Antenna):
    """
    Implements the antenna pattern for the Radio Astronomy Service 
    according to recommendation ITU-R SA.509-3.
    """
    
    def __init__(self, param: ParametersRas):
        super().__init__()
        # Set basic attributes
        self.diameter = param.diameter
        self.efficiency = param.antenna_efficiency
        self.wavelength = param.SPEED_OF_LIGHT/(param.frequency*1e6)
        
        # Diagram parameters
        self.g_0 = 10*np.log10(self.efficiency*\
                               (np.pi*self.diameter/self.wavelength)**2)
        self.phi_0 = 20*np.sqrt(3)/(self.diameter/self.wavelength)
        
        # Limit parameters
        self.phi_1 = self.phi_0*np.sqrt(20/3)
        self.phi_2 = 10**((49-self.g_0)/25)
        
    def calculate_gain(self, *args, **kwargs) -> np.array:
        phi = np.absolute(kwargs["phi_vec"])
        
        gain = np.zeros_like(phi)
        
        # First part
        interval_idx = np.where(np.logical_and(phi >= 0, phi < self.phi_1))
        gain[interval_idx] = self.g_0 - 3*(phi[interval_idx]/self.phi_0)**2
        # Second part
        interval_idx = np.where(np.logical_and(phi >= self.phi_1, phi < self.phi_2))
        gain[interval_idx] = self.g_0 - 20
        # Third part
        interval_idx = np.where(np.logical_and(phi >= self.phi_2, phi < 48))
        gain[interval_idx] = 29 - 25*np.log10(phi(interval_idx))
        # Fourth part
        interval_idx = np.where(np.logical_and(phi >= 48, phi < 80))
        gain[interval_idx] = -13
        # Fifth part
        interval_idx = np.where(np.logical_and(phi >= 80, phi < 120))
        gain[interval_idx] = -8
        # Sixth part
        interval_idx = np.where(np.logical_and(phi >= 120, phi < 180))
        gain[interval_idx] = -13
        
if __name__ == '__main__':
    import matplotlib.pyplot as plt

    phi = np.linspace(0.1, 180, num = 100000)