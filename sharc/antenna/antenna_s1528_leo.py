# -*- coding: utf-8 -*-
"""
Created on Fri Marc 12 15:58:01 2021

@author: Luciano Camilo
"""
from sharc.antenna.antenna import Antenna
from sharc.parameters.parameters_fss_ss import ParametersFssSs

import numpy as np


class AntennaS1528_LEO(Antenna):
    """
    Implements Recommendation ITU-R S.1528-0: Satellite antenna radiation
    patterns for LEO orbit satellite antennas operating in the
    fixed-satellite service below 30 GHz
    """

    def __init__(self, param: ParametersFssSs):
        super().__init__()
        self.peak_gain = param.antenna_gain

        # near-in-side-lobe level (dB) relative to the peak gain required by
        # the system design
        self.l_s = param.antenna_l_s
        # for elliptical antennas, this is the ratio major axis/minor axis
        # we assume circular antennas, so z = 1
        self.z = 1
        # far-out side-lobe level [dBi]
        self.l_f = 5
        self.psi_b = param.antenna_3_dB*2
        self.b = 1.6
        self.alpha = 2
        self.y = 1.5*self.b
        self.z = 20.4

    def calculate_gain(self, *args, **kwargs) -> np.array:
        psi = np.absolute(kwargs["off_axis_angle_vec"])
        gain = np.zeros(len(psi))

        idx_0 = np.where((psi <= self.b))[0]
        gain[idx_0] = self.peak_gain

        idx_1 = np.where(((self.b) < psi) & (psi <= self.y))[0]
        gain[idx_1] = self.peak_gain - 3 * np.power(psi[idx_1] / self.psi_b, self.alpha)

        idx_2 = np.where((self.y < psi) & (psi <= self.z))[0]
        gain[idx_2] = self.peak_gain + self.l_s - 25 * np.log10(psi[idx_2]/self.y)

        idx_3 = np.where((self.z < psi) & (psi <= 180))[0]
        gain[idx_3] = self.l_f

        return gain


if __name__ == '__main__':
    import matplotlib.pyplot as plt

    # initialize antenna parameters
    param = ParametersFssSs()
    param.antenna_gain = 35
    param.antenna_pattern = "ITU-R S.1528-LEO"
    param.antenna_3_dB = 1.6
    psi = np.linspace(0, 90, num=10000)

    param.antenna_l_s = -6.75
    antenna = AntennaS1528_LEO(param)
    gain12 = antenna.calculate_gain(off_axis_angle_vec=psi)

    fig = plt.figure(figsize=(8, 7), facecolor='w', edgecolor='k')  # create a figure object
    plt.plot(psi, gain12, "-b", label="$L_S = -6.75$ dB")

    # plt.ylim((-40, 10))
    plt.xlim((0, 90))
    plt.title("ITU-R S.1528-0 LEO antenna radiation pattern")
    plt.xlabel("Relative off-axis angle, $\psi/\psi_{3dB}$")
    plt.ylabel("Gain relative to $G_{max}$ [dB]")
    plt.legend(loc="upper right")

#    ax = plt.gca()
#    ax.set_yticks([-30, -20, -10, 0])
#    ax.set_xticks(np.linspace(1, 9, 9).tolist() + np.linspace(10, 100, 10).tolist())

    plt.grid()
    plt.show()
