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
        
    def calculate_gain(self):
        pass