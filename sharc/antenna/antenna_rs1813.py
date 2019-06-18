# -*- coding: utf-8 -*-
"""
Created on Thu Feb 21 11:22:08 2019

@author: Edgar Souza
"""

from sharc.antenna.antenna import Antenna
from sharc.parameters.parameters_eess_passive import ParametersEessPassive

import numpy as np
import math

class AntennaRS1813(Antenna):
    """
    Implements the reference antenna pattern described in 
    Recommendation ITU-R RS.1813-1. This is the antenna pattern for 
    Earth exploration-satellite service (EESS) passive sensors to be 
    used in compatibility studies in the frequency range 1.4-100 GHz.
    
    This implementation is in accordance with recommends 2, which refers to
    the case where a few interference sources dominate, or where peak 
    interference values are required in the analysis, the following equations 
    for the antenna pattern for spaceborne passive sensors should be used, 
    for antenna diameters greater than 2 times the wavelength.
    """

    def __init__(self, param: ParametersEessPassive):
        super().__init__()
        self.lmbda = 3e8 / ( param.frequency * 1e6 )
        self.d_lmbda = param.antenna_diameter / self.lmbda
        
        # for sensor F3 with n = 60% and D = 2.2 m, G_max = 52.7 dBi
        self.peak_gain = 10 * math.log10(param.antenna_efficiency * math.pow(math.pi * self.d_lmbda, 2))
        #self.peak_gain = param.antenna_gain
        
        self.phi_m = 22 / self.d_lmbda * \
                math.sqrt(5.5 + 5 * math.log10(math.pow(param.antenna_efficiency, 2) * self.d_lmbda))


    def calculate_gain(self, *args, **kwargs) -> np.array:
        phi = np.absolute(kwargs["off_axis_angle_vec"])

        gain = np.zeros(phi.shape)

        id0 = np.where(phi <= 69)[0]
        gain[id0] = self.peak_gain - 0.0018 * np.power(self.d_lmbda * phi[id0], 2)
        
        id1 = np.where((self.phi_m < phi) & (phi <= 69))[0]
        gain[id1] = np.maximum( gain[id1], 33 - 5*math.log10(self.d_lmbda) - 25*np.log10(phi[id1]) )
        
        id2 = np.where((69 < phi) & (phi <= 180))[0]
        gain[id2] = -13 - 5*math.log10(self.d_lmbda)
        
        gain = np.maximum( -23, gain )

        return gain


if __name__ == '__main__':
    import matplotlib.pyplot as plt

    phi = np.linspace(0.1, 180, num = 10000)

    # initialize antenna parameters
    param = ParametersEessPassive()
    param.antenna_pattern = "ITU-R RS.1813-1"
    param.frequency = 23900
    param.antenna_gain = 52
    param.antenna_diameter = 2.2
    param.antenna_efficiency = 0.6
    antenna = AntennaRS1813(param)

    gain = antenna.calculate_gain(off_axis_angle_vec = phi)

    fig = plt.figure(figsize=(8,7), facecolor='w', edgecolor='k')  # create a figure object
    plt.semilogx(phi, gain - param.antenna_gain, "-b", label = "$f = 23.9$ GHz")

    plt.title("ITU-R RS.1813-1 antenna radiation pattern")
    plt.xlabel("Off-axis angle $\phi$ [deg]")
    plt.ylabel("Normalized antenna gain [dBi]")
    plt.legend(loc="lower left")
    plt.xlim((phi[0], phi[-1]))
    plt.ylim((-80, 10))

    plt.grid()
    plt.show()
