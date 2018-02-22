# -*- coding: utf-8 -*-
"""
Created on Thu Aug 24 18:41:17 2017

@author: edgar
"""

from sharc.antenna.antenna import Antenna
from sharc.parameters.parameters_fss_es import ParametersFssEs

import numpy as np

class AntennaS580(Antenna):
    """
    Implements the Earth station antenna pattern in the EESS/ISS service
    according to Recommendation ITU-R S.580-6
    """

    def __init__(self, param: ParametersFssEs):
        super().__init__()
        self.peak_gain = param.antenna_gain
        lmbda = 3e8 / ( param.frequency * 1e6 )

        self.phi_min = 1
        if 100 * lmbda / param.diameter > 1:
           self.phi_min = 100 * lmbda / param.diameter

    def calculate_gain(self, *args, **kwargs) -> np.array:
        phi = np.absolute(kwargs["off_axis_angle_vec"])

        gain = np.zeros(phi.shape)

        idx_0 = np.where(phi < self.phi_min)[0]
        gain[idx_0] = self.peak_gain

        idx_1 = np.where((self.phi_min <= phi) & (phi < 20))[0]
        gain[idx_1] = 29 - 25 * np.log10(phi[idx_1])

        idx_2 = np.where((20 <= phi) & (phi <= 180))[0]
        gain[idx_2] = -10

        return gain


if __name__ == '__main__':
    import matplotlib.pyplot as plt

    phi = np.linspace(0.1, 100, num = 100000)

    # initialize antenna parameters
    param27 = ParametersFssEs()
    param27.antenna_pattern = "ITU-R S.580-6"
    param27.frequency = 27000
    param27.antenna_gain = 50
    param27.diameter = 9.6
    antenna27 = AntennaS580(param27)

    gain27 = antenna27.calculate_gain(phi_vec=phi)

    param = ParametersFssEs()
    param.antenna_pattern = "ITU-R S.580-6"
    param.frequency = 27000
    param.antenna_gain = 50
    param.diameter = 0.45
    antenna = AntennaS580(param)
    gain = antenna.calculate_gain(phi_vec=phi)

    fig = plt.figure(figsize=(8,7), facecolor='w', edgecolor='k')  # create a figure object

    plt.semilogx(phi, gain27 - param27.antenna_gain, "-b", label = "$f = 27$ $GHz,$ $D = 9.6$ $m$")
    plt.semilogx(phi, gain - param.antenna_gain, "-r", label = "$f = 27$ $GHz,$ $D = 0.45$ $m$")

    plt.title("ITU-R S.580 antenna radiation pattern")
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
