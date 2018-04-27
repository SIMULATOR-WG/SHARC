# -*- coding: utf-8 -*-
"""
Created on Wed Apr  4 17:08:00 2018
@author: Calil
"""

from sharc.antenna.antenna import Antenna
from sharc.parameters.parameters_fss_es import ParametersFssEs

import numpy as np
import math

class AntennaF1245(Antenna):
    """
    Implements reference radiation patterns for HAPS gateway antennas
    for use in coordination studies and interference assessment in the
    frequency range from 1 GHz to about 70 GHz. (ITU-R F.1245-2)
    """

    def __init__(self, param: ParametersFssEs):
        super().__init__()
        self.peak_gain = param.antenna_gain
        lmbda = 299792458 / ( param.frequency * 1e6 )
        self.d_lmbda = param.diameter / lmbda

        self.g_l = 2 + 15 * math.log10(self.d_lmbda)
        self.phi_m = (20 / self.d_lmbda) * math.sqrt(self.peak_gain - self.g_l)
        self.phi_r = 12.02 * math.pow(self.d_lmbda, -0.6)

    def calculate_gain(self, *args, **kwargs) -> np.array:
        off_axis = np.absolute(kwargs["off_axis_angle_vec"])
        
        if self.d_lmbda > 100:
            gain = self.calculate_gain_greater(off_axis)
        else:
            gain = self.calculate_gain_less(off_axis)
            
        return gain

    def calculate_gain_greater(self, phi: float) -> np.array:
        """
        For frequencies in the range 1 GHz to about 70 GHz, in cases where the
        ratio between the antenna diameter and the wavelength is GREATER than
        100, this method should be used.
        Parameter
        ---------
            phi : off-axis angle [deg]
        Returns
        -------
            a numpy array containing the gains in the given angles
        """
        gain = np.zeros(phi.shape)

        idx_0 = np.where(phi < self.phi_m)[0]
        gain[idx_0] = self.peak_gain - 2.5e-3 * np.power(self.d_lmbda * phi[idx_0], 2)

        phi_thresh = max(self.phi_m,self.phi_r)
        idx_1 = np.where((self.phi_m <= phi) & (phi < phi_thresh))[0]
        gain[idx_1] = self.g_l

        idx_2 = np.where((phi_thresh <= phi) & (phi < 48))[0]
        gain[idx_2] = 29 - 25 * np.log10(phi[idx_2])

        idx_3 = np.where((48 <= phi) & (phi <= 180))[0]
        gain[idx_3] = -13

        return gain

    def calculate_gain_less(self, phi: float) -> np.array:
        """
        For frequencies in the range 1 GHz to about 70 GHz, in cases where the
        ratio between the antenna diameter and the wavelength is LESS than
        or equal to 100, this method should be used.
        Parameter
        ---------
            phi : off-axis angle [deg]
        Returns
        -------
            a numpy array containing the gains in the given angles
        """
        gain = np.zeros(phi.shape)

        idx_0 = np.where(phi < self.phi_m)[0]
        gain[idx_0] = self.peak_gain - 0.0025 * np.power(self.d_lmbda * phi[idx_0], 2)

        idx_1 = np.where((self.phi_m <= phi) & (phi < 48))[0]
        gain[idx_1] = 39 - 5 * math.log10(self.d_lmbda) - 25 * np.log10(phi[idx_1])

        idx_2 = np.where((48 <= phi) & (phi < 180))[0]
        gain[idx_2] = -3 - 5 * math.log10(self.d_lmbda)

        return gain
    

if __name__ == '__main__':
    import matplotlib.pyplot as plt

    off_axis_angle = np.linspace(0.1, 180, num = 100000)

    # initialize antenna parameters
    param = ParametersFssEs()
    param.frequency = 27250
    param.antenna_gain = 53.3
    param.diameter = 2
    antenna = AntennaF1245(param)

    gain = antenna.calculate_gain(off_axis_angle_vec = off_axis_angle)
    
    fig = plt.figure(figsize=(6,5), facecolor='w', edgecolor='k')  # create a figure object

    plt.semilogx(off_axis_angle, gain, "-b")

    #plt.title("ITU-R F.1245-2 antenna radiation pattern")
    plt.xlabel("Off-axis angle $\phi$ [deg]")
    plt.ylabel("Gain [dBi]")
    plt.xlim((off_axis_angle[0], off_axis_angle[-1]))
    plt.ylim((-20, 60))

    plt.grid()
    plt.show()