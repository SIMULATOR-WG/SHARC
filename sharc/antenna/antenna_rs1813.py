#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 16:45:35 2017

@author: carlosrodriguez
"""

from sharc.antenna.antenna import Antenna
from sharc.parameters.parameters_fss_ss import ParametersFssSs

import numpy as np

class AntennaRS1813(Antenna):
    """
    Implements the Earth station antenna pattern in the EESS/ISS service
    according to Recommendation ITU-R RS.1813-1
    Reference antenna pattern for passive sensors operating in the 
    Earth exploration-satellite service (passive) to be used in 
    compatibility analyses in the frequency range 1.4-100 GHz
    
    Implementation of item 1) of recommendation 
    """
    
    def __init__(self, param: ParametersFssSs):
        super().__init__()
        self.lmbda = 3e8 / ( param.frequency * 1e6 )
        self.peak_gain = 10 * np.log10(60 * np.power(3.141516,2) * np.power(param.diameter / self.lmbda,2))

        self.phi_min = 1
        #efficiency assumed as 60%
        self.phi_min = 22 * self.lmbda / param.diameter * np.sqrt(5.5 + 5 * np.log10(60 * 60 * param.diameter / self.lmbda))
        self.d_to_lmbda = param.diameter / self.lmbda
        
    def calculate_gain(self, *args, **kwargs) -> np.array:
        phi = np.absolute(kwargs["phi_vec"])
        
        gain = np.zeros(phi.shape)
                
        idx_0 = np.where(phi < self.phi_min)[0]
        gain[idx_0] = self.peak_gain - 0.0018 * np.power(phi[idx_0] * self.d_to_lmbda, 2)
        
        idx_1 = np.where((self.phi_min <= phi) & (phi < 69))[0]
        gain[idx_1] = np.maximum([self.peak_gain - 0.0018 * np.power(phi[idx_1] * self.d_to_lmbda, 2)],[ 33 - 5 * np.log10(self.d_to_lmbda) - 25 * np.log10(phi[idx_1])])
            
        idx_2 = np.where((69 <= phi) & (phi <= 180))[0]
        gain[idx_2] = -13 - 5 * np.log10(self.d_to_lmbda)
        
        return gain
        
        
if __name__ == '__main__':
    import matplotlib.pyplot as plt

    phi = np.linspace(0.1, 10, num = 100000)
    
    # initialize antenna parameters
    param26 = ParametersFssSs()
    param26.antenna_pattern = "ITU-R RS.1813-1"
    param26.frequency = 23800
    param26.diameter = 0.6
    antenna26 = AntennaRS1813(param26)

    gain26 = antenna26.calculate_gain(phi_vec=phi)
    
    fig = plt.figure(figsize=(8,7), facecolor='w', edgecolor='k')  # create a figure object
    
    plt.plot(phi, gain26 - 61.276437077350174, "-b", label = "$f = 23.8$ $GHz,$ $D = 0.6$ $m$")
    
    plt.title("ITU-R RS.1813-1 antenna radiation pattern")
    plt.xlabel("Off-axis angle $\phi$ [deg]")
    plt.ylabel("Gain relative to $G_m$ [dB]")
    plt.legend(loc="lower left")
    plt.xlim((phi[0], phi[-1]))
    plt.ylim((-100, 10))
    
    #ax = plt.gca()
    #ax.set_yticks([-30, -20, -10, 0])
    #ax.set_xticks(np.linspace(1, 9, 9).tolist() + np.linspace(10, 100, 10).tolist())

    plt.grid()
    plt.show()        