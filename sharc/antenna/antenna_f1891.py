# -*- coding: utf-8 -*-
"""
Created on Thu Oct 19 12:24:13 2017

@author: edgar
"""

from sharc.antenna.antenna import Antenna
from sharc.parameters.parameters_haps import ParametersHaps

import numpy as np

class AntennaF1891(Antenna):
    """
    Implements the antenna radiation pattern described in ITU-R F.1891. This is
    a phased array antenna that is described in, and that complies with,
    Resolution 221 (Rev.WRC-07) and can be used in both the HAPS gateway
    (ground) station and in the HAPS (airborne) platform.
    """

    def __init__(self, param: ParametersHaps):
        super().__init__()
        self.peak_gain = param.antenna_gain
        self.psi_b = np.sqrt(7442/np.power(10, 0.1 * self.peak_gain))
        self.l_n = param.antenna_l_n
        self.l_f = self.peak_gain - 73
        self.psi_1 = self.psi_b * np.sqrt(-self.l_n/3)
        self.psi_2 = 3.745 * self.psi_b
        self.x = self.peak_gain + self.l_n + 60*np.log10(self.psi_2)
        self.psi_3 = np.power(10, (self.x - self.l_f)/60)

    def calculate_gain(self, *args, **kwargs) -> np.array:
        psi = np.absolute(kwargs["off_axis_angle_vec"])

        gain = self.l_f * np.ones(len(psi))

        idx_0 = np.where(psi <= self.psi_1)[0]
        gain[idx_0] = self.peak_gain - 3*np.power(psi[idx_0]/self.psi_b, 2)

        idx_1 = np.where((self.psi_1 < psi) & (psi <= self.psi_2))[0]
        gain[idx_1] = self.peak_gain + self.l_n

        idx_2 = np.where((self.psi_2 < psi) & (psi <= self.psi_3))[0]
        gain[idx_2] = self.x - 60*np.log10(psi[idx_2])

        return gain


if __name__ == '__main__':
    import matplotlib.pyplot as plt

    # initialize antenna parameters
    param = ParametersHaps()
    param.antenna_pattern = "ITU-R F.1891"
    param.antenna_l_n = -25
    psi = np.linspace(0, 90, num = 10001)

    param.antenna_gain = 47
    antenna47 = AntennaF1891(param)
    gain47 = antenna47.calculate_gain(phi_vec=psi)

    param.antenna_gain = 30
    antenna30 = AntennaF1891(param)
    gain30 = antenna30.calculate_gain(phi_vec=psi)

    fig = plt.figure(figsize=(16,6), facecolor='w', edgecolor='k')  # create a figure object
    ax1 = fig.add_subplot(121)

    ax1.plot(psi, gain47, "-b", label="$G_m = 47$ dB")
    ax1.set_ylim((-40, 60))
    ax1.set_xlim((0, 20))
    ax1.set_title("ITU-R F.1891 antenna radiation pattern")
    ax1.set_xlabel("Off-axis angle [deg]")
    ax1.set_ylabel("Gain [dBi]")
    ax1.legend(loc="upper right")
    ax1.set_yticks([-40, -20, 0, 20, 40, 60])
    ax1.set_xticks(np.linspace(0, 20, 11).tolist())
    ax1.grid(True)

    ax2 = fig.add_subplot(122)
    ax2.plot(psi, gain30, "-r", label="$G_m = 30$ dB")
    ax2.set_ylim((-60, 40))
    ax2.set_xlim((0, 90))
    ax2.set_title("ITU-R F.1891 antenna radiation pattern")
    ax2.set_xlabel("Off-axis angle [deg]")
    ax2.set_ylabel("Gain [dBi]")
    ax2.legend(loc="upper right")
    ax2.set_yticks([-60, -40, -20, 0, 20, 40])
    ax2.set_xticks(np.linspace(0, 90, 10).tolist())
    ax2.grid(True)

    plt.show()
