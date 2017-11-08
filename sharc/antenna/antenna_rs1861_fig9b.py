#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  8 13:40:27 2017

@author: carlosrodriguez
"""

from sharc.antenna.antenna import Antenna
from sharc.parameters.parameters_fss_es import ParametersFssEs

import numpy as np

class AntennaRS1861FIG9b(Antenna):
    """
    Implements the EESS sensor antenna for 23.6-24 GHz for Sensor type F2
    """
    
    def __init__(self, param: ParametersFssEs):
        super().__init__()
        self.peak_gain = param.antenna_gain
        

        self.phi_1 = 1.5       
        self.phi_2 = 10
                
    def calculate_gain(self, *args, **kwargs) -> np.array:
        phi = np.absolute(kwargs["phi_vec"])
        
        gain = np.zeros(phi.shape)
        "Form parabolic equation according to figure 9b, a*attˆ2+b*att+angle = 0"
        a = 0.00055644
        b = 0.06112878
        
        "From line eq for second part of fig 9b att=m*angle+c"
        m = -3.2944176
        c = -32.058824

        idx_0 = np.where((0 <= phi) & (phi < self.phi_1))[0]
        gain[idx_0] = self.peak_gain + ( - 1 * b + np.sqrt(np.power(b,2) - 4 * a * phi[idx_0]))/(2 * a)

        idx_1 = np.where((self.phi_1 <= phi) & (phi <= self.phi_2))[0]
        gain[idx_1] = self.peak_gain + m * phi[idx_1] + c
        
        idx_2 = np.where((self.phi_2 <= phi) & (phi < 18))[0]
        gain[idx_2] = self.peak_gain - 65
         
        return gain
        
        
if __name__ == '__main__':
    import matplotlib.pyplot as plt

    phi = np.linspace(0.1, 100, num = 100000)
    
    # initialize antenna parameters
    param27 = ParametersFssEs()
    param27.antenna_pattern = "ITU-R RS.1861 Figure 9b"
    param27.frequency = 26000
    param27.antenna_gain = 50
    antenna27 = AntennaRS1861FIG9b(param27)

    gain27 = antenna27.calculate_gain(phi_vec=phi)
    

    fig = plt.figure(figsize=(8,7), facecolor='w', edgecolor='k')  # create a figure object
    
    plt.plot(phi, gain27 - param27.antenna_gain, "-b", label = "$f = 26$ $GHz$")

    plt.title("ITU-R RS.1861 Figure 9b antenna radiation pattern")
    plt.xlabel("Off-axis angle $\phi$ [deg]")
    plt.ylabel("Gain relative to $G_m$ [dB]")
    plt.legend(loc="lower left")
    plt.xlim((phi[0], 18))
    plt.ylim((-80, 10))
    
    #ax = plt.gca()
    #ax.set_yticks([-30, -20, -10, 0])
    #ax.set_xticks(np.linspace(1, 9, 9).tolist() + np.linspace(10, 100, 10).tolist())

    plt.grid()
    plt.show()        