# -*- coding: utf-8 -*-
"""
Created on Tue Dec 15 11:53:58 2020

@modified: Luciano Camilo Tue Jan 26 13:49:25 2021
"""

import numpy as np
from sharc.parameters.parameters_fs import ParametersFs


class AntennaOmniF1336(object):
    """
    Implements a omni-directional antenna pattern of an FS antenna array following ITU-R F.1336-5, Annex I

    Attributes
    ----------
        antenna_gain (float): maximum gain of FS omni antenna
    """

    def __init__(self, param: ParametersFs):
        """
        Constructs an AntennaElementImt object.

        Parameters
        ---------
            param (ParametersFS): antenna FS parameters
        """
        super().__init__()
        self.g_max = param.antenna_gain
        # 3 dB beamwidth calculation
        self.theta_3db = 1 / ((((10 ** (self.g_max / 10)) + 172.4) / 191) ** 2 - 0.818)
        self.m = -3 / (10 * np.log10(np.cos(np.radians((1 / 2) * self.theta_3db))))
        # Antenna number of elements (for comparison with practical antenna models)
        self.n = 2.7832 / ((2 * np.pi * (3 / 4)) * (np.sin(np.radians(self.theta_3db / 2))))

    def calculate_gain(self, **kwargs) -> np.array:

        theta = np.asarray(kwargs["theta_vec"])
        g1 = self.g_max - 12 * ((theta / self.theta_3db) ** 2)
        # k - increase side-lobes above what would be expected for an antenna with improved side-lobe performance
        k = 0
        g2 = self.g_max - 12 + 10 * np.log10((np.maximum((abs(theta) / self.theta_3db), 1) ** -1.5) + k)
        gain = np.maximum(g1, g2)
        # For normalized pattern test please put the maximum gain in Gn variable and uncomment the code
        # Gn=13
        # gain = gain - Gn
        # print(gain)
        return gain


if __name__ == '__main__':

    from sharc.parameters.parameters_fs import ParametersFs
    from matplotlib import pyplot as plt
    """
    Test routine - Comparison between SHARC code and Annex I of ITU R - F.1336-5
    """

    param = ParametersFs()
    param.antenna_gain = 13
    antenna = AntennaOmniF1336(param)
    phi = np.arange(0, 360, step=1)
    theta = np.arange(0, 90, step=1)
    # theta = np.arange(-90, 90, step=1)

    plt.figure(1)
    plt.legend(loc='upper left')
    plt.plot(theta, antenna.calculate_gain(off_axis_angle_vec=phi, theta_vec=theta), 'b--',
             label='ITU-R F-1336-5 Omni Pattern')

    plt.xlim(0, 90)
    # plt.xlim(-90, 90)
    plt.grid(which='minor', alpha=0.2)
    plt.grid(which='major', alpha=0.5)
    plt.grid(True, color='b', linestyle='--', linewidth=0.2)
    plt.title('ITU-R F-1336-5 Omni Pattern')
    plt.xlabel('Elevation angle (degrees)')
    plt.ylabel('Gain (dBi)')
    plt.legend()
    plt.show()
