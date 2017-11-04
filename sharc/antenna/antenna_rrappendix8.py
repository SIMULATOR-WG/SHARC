#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov  4 15:28:30 2017

@author: carlosrodriguez
"""


from sharc.antenna.antenna import Antenna
from sharc.parameters.parameters_fss_es import ParametersFssEs

import numpy as np

class AntennaRRAppendix8(Antenna):
    """
    Implements the Earth station antenna pattern in the EESS service
    according to RR Appendix 8 for GSO networks and D/lambda >100
    """
    
    def __init__(self, param: ParametersFssEs):
        super().__init__()
        self.peak_gain = param.antenna_gain
        self.Dlmbda = np.power(10, (self.peak_gain - 7.7) / 20)
        self.G1 = 2 + 15 * np.log10(self.Dlmbda)

        self.phi_m = 1        
        self.phi_r = 1
        self.phi_m = 20 * 1 / self.Dlmbda * np.sqrt(self.peak_gain - self.G1) 
        self.phi_r = 15.85 * np.power(self.Dlmbda,-0.6) 
        
    def calculate_gain(self, *args, **kwargs) -> np.array:
        phi = np.absolute(kwargs["phi_vec"])
        
        gain = np.zeros(phi.shape)
        
        if self.Dlmbda >= 100:
            idx_0 = np.where((0 <= phi) & (phi < self.phi_m))[0]
            gain[idx_0] = self.peak_gain - 0.0025 * np.power(self.Dlmbda * phi[idx_0], 2)

            idx_1 = np.where((self.phi_m <= phi) & (phi <= self.phi_r))[0]
            gain[idx_1] = self.G1
        
            idx_2 = np.where((self.phi_r <= phi) & (phi < 48))[0]
            gain[idx_2] = 32 - 25 * np.log10(phi[idx_2])

            idx_3 = np.where((48 <= phi) & (phi <= 180))[0]
            gain[idx_3] = -10
            
        elif self.Dlmbda < 100 :
            idx_0 = np.where((0 <= phi) & (phi < self.phi_m))[0]
            gain[idx_0] = self.peak_gain - 0.0025 * np.power(self.Dlmbda * phi[idx_0], 2)

            idx_1 = np.where((self.phi_m <= phi) & (phi <= self.phi_r))[0]
            gain[idx_1] = self.G1
        
            idx_2 = np.where((self.phi_r <= phi) & (phi < 48))[0]
            gain[idx_2] = 52 - 10 * np.log10(self.Dlmbda) - 25 * np.log10(phi[idx_2])

            idx_3 = np.where((48 <= phi) & (phi <= 180))[0]
            gain[idx_3] = -10 - 10 * np.log10(self.Dlmbda)
        return gain
        
        
if __name__ == '__main__':
    import matplotlib.pyplot as plt

    phi = np.linspace(0.1, 100, num = 100000)
    
    # initialize antenna parameters
    param27 = ParametersFssEs()
    param27.antenna_pattern = "RR Appendix 8"
    param27.frequency = 26000
    param27.antenna_gain = 50
    antenna27 = AntennaRRAppendix8(param27)

    gain27 = antenna27.calculate_gain(phi_vec=phi)
    

    fig = plt.figure(figsize=(8,7), facecolor='w', edgecolor='k')  # create a figure object
    
    plt.semilogx(phi, gain27 - param27.antenna_gain, "-b", label = "$f = 26$ $GHz$")

    plt.title("RR Appendix 8 antenna radiation pattern")
    plt.xlabel("Off-axis angle $\phi$ [deg]")
    plt.ylabel("Gain relative to $G_m$ [dB]")
    plt.legend(loc="lower left")
    plt.xlim((phi[0], phi[-1]))
    plt.ylim((-80, 10))
    
    #ax = plt.gca()
    #ax.set_yticks([-30, -20, -10, 0])
    #ax.set_xticks(np.linspace(1, 9, 9).tolist() + np.linspace(10, 100, 10).tolist())

    plt.grid()
    plt.show()        