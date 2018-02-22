# -*- coding: utf-8 -*-
"""
Created on Thu Nov 16 17:34:11 2017

@author: Calil
"""

from sharc.antenna.antenna import Antenna
from sharc.parameters.parameters_ras import ParametersRas

import numpy as np

class AntennaSA509(Antenna):
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

        # Effective area
        self.effective_area = self.efficiency*(np.pi*self.diameter**2)/4

        # Diagram parameters
        self.g_0 = 10*np.log10(self.efficiency*\
                               (np.pi*self.diameter/self.wavelength)**2)
        self.phi_0 = 20*np.sqrt(3)/(self.diameter/self.wavelength)

        # Limit parameters
        self.phi_1 = self.phi_0*np.sqrt(20/3)
        self.phi_2 = 10**((49-self.g_0)/25)

    def calculate_gain(self, *args, **kwargs) -> np.array:
        phi = np.absolute(kwargs["off_axis_angle_vec"])

        gain = np.zeros_like(phi)

        # First part
        interval_idx = np.where(np.logical_and(phi >= 0, phi < self.phi_1))
        gain[interval_idx] = self.g_0 - 3*(phi[interval_idx]/self.phi_0)**2
        # Second part
        interval_idx = np.where(np.logical_and(phi >= self.phi_1, phi < self.phi_2))
        gain[interval_idx] = self.g_0 - 20
        # Third part
        interval_idx = np.where(np.logical_and(phi >= self.phi_2, phi < 48))
        gain[interval_idx] = 29 - 25*np.log10(phi[interval_idx])
        # Fourth part
        interval_idx = np.where(np.logical_and(phi >= 48, phi < 80))
        gain[interval_idx] = -13
        # Fifth part
        interval_idx = np.where(np.logical_and(phi >= 80, phi < 120))
        gain[interval_idx] = -8
        # Sixth part
        interval_idx = np.where(np.logical_and(phi >= 120, phi <= 180))
        gain[interval_idx] = -13

        return gain

if __name__ == '__main__':
    import matplotlib.pyplot as plt

    par = ParametersRas();
    par.diameter = 1
    par.antenna_efficiency = 1
    par.frequency = 43000
    par.SPEED_OF_LIGHT = 3e8

    antenna1 = AntennaSA509(par)
    par.diameter = 7
    antenna7 = AntennaSA509(par)
    par.diameter = 10
    antenna10 = AntennaSA509(par)

    phi = np.linspace(0.1, 180, num = 100000)
    gain1  = antenna1.calculate_gain(phi_vec = phi)
    gain7  = antenna7.calculate_gain(phi_vec = phi)
    gain10  = antenna10.calculate_gain(phi_vec = phi)

    fig = plt.figure(figsize=(8,7), facecolor='w', edgecolor='k')  # create a figure object

    plt.semilogx(phi, gain1, "-b", label = "$f = 43$ $GHz,$ $D = 1$ $m$")
    plt.semilogx(phi, gain7, "-r", label = "$f = 43$ $GHz,$ $D = 7$ $m$")
    plt.semilogx(phi, gain10, "-k", label = "$f = 43$ $GHz,$ $D = 10$ $m$")

    plt.title("ITU-R SA.509-3 antenna radiation pattern")
    plt.xlabel("Off-axis angle $\phi$ [deg]")
    plt.ylabel("Gain [dBi]")
    plt.legend(loc="lower left")
    plt.xlim((phi[0], phi[-1]))

    plt.grid()
    plt.show()

