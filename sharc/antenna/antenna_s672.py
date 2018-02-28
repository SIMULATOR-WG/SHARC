# -*- coding: utf-8 -*-
"""
Created on Thu Apr 13 17:18:59 2017

@author: edgar
"""

from sharc.antenna.antenna import Antenna
from sharc.parameters.parameters_fss_ss import ParametersFssSs

import numpy as np
import sys

class AntennaS672(Antenna):
    """
    Implements the satellite antenna pattern in the fixed-satellite service
    according to Recommendation ITU-R S.672-4 Annex 1
    """

    def __init__(self, param: ParametersFssSs):
        super().__init__()
        self.peak_gain = param.antenna_gain
        self.l_s = param.antenna_l_s

        if self.l_s == -20:
            self.a = 2.58
        elif self.l_s == -25:
            self.a = 2.88
        elif self.l_s == -30:
            self.a = 3.16
        else:
            sys.stderr.write("ERROR\nInvalid AntennaS672 L_s parameter: " + self.l_s)
            sys.exit(1)

        self.b = 6.32

        self.psi_0 = param.antenna_3_dB/2
        self.psi_1 = self.psi_0 * np.power(10, (self.peak_gain + self.l_s + 20)/25)


    def calculate_gain(self, *args, **kwargs) -> np.array:
        psi = np.absolute(kwargs["off_axis_angle_vec"])

        gain = np.zeros(len(psi))

        idx_0 = np.where(psi < self.psi_0)
        gain[idx_0] = self.peak_gain

        idx_1 = np.where((self.psi_0 <= psi) & (psi <= self.a * self.psi_0))[0]
        gain[idx_1] = self.peak_gain - 3 * np.power(psi[idx_1]/self.psi_0, 2)

        idx_2 = np.where((self.a * self.psi_0 < psi) & (psi <= self.b * self.psi_0))[0]
        gain[idx_2] = self.peak_gain + self.l_s

        idx_3 = np.where((self.b * self.psi_0 < psi) & (psi <= self.psi_1))[0]
        gain[idx_3] = self.peak_gain + self.l_s + 20 - 25 * np.log10(psi[idx_3]/self.psi_0)

        return gain


if __name__ == '__main__':
    import matplotlib.pyplot as plt

    # initialize antenna parameters
    param = ParametersFssSs()
    param.antenna_gain = 50
    param.antenna_pattern = "ITU-R S.672-4"
    param.antenna_3_dB = 2
    psi = np.linspace(1, 30, num = 1000)

    param.antenna_l_s = -20
    antenna = AntennaS672(param)
    gain20 = antenna.calculate_gain(phi_vec=psi)

    param.antenna_l_s = -25
    antenna = AntennaS672(param)
    gain25 = antenna.calculate_gain(phi_vec=psi)

    param.antenna_l_s = -30
    antenna = AntennaS672(param)
    gain30 = antenna.calculate_gain(phi_vec=psi)

    fig = plt.figure(figsize=(12,7), facecolor='w', edgecolor='k')  # create a figure object

    plt.semilogx(psi, gain20 - param.antenna_gain, "-b", label="$L_S = -20$ dB")
    plt.semilogx(psi, gain25 - param.antenna_gain, "-r", label="$L_S = -25$ dB")
    plt.semilogx(psi, gain30 - param.antenna_gain, "-g", label="$L_S = -30$ dB")

    plt.ylim((-33.8, 0))
    plt.xlim((1, 100))
    plt.title("ITU-R S.672-4 antenna radiation pattern")
    plt.xlabel("Relative off-axis angle, $\psi/\psi_0$")
    plt.ylabel("Gain relative to $G_m$ [dB]")
    plt.legend(loc="upper right")

    ax = plt.gca()
    ax.set_yticks([-30, -20, -10, 0])
    ax.set_xticks(np.linspace(1, 9, 9).tolist() + np.linspace(10, 100, 10).tolist())

    plt.grid()
    plt.show()
