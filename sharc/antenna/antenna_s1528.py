# -*- coding: utf-8 -*-
"""
Created on Tue Aug 15 14:49:01 2017

@author: edgar
"""

from sharc.antenna.antenna import Antenna
from sharc.parameters.parameters_fss_ss import ParametersFssSs

import math
import numpy as np
import sys

class AntennaS1528(Antenna):
    """
    Implements Recommendation ITU-R S.1528-0: Satellite antenna radiation
    patterns for non-geostationary orbit satellite antennas operating in the
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
        self.l_f = 0

        # back-lobe level
        self.l_b = np.maximum(0, 15 + self.l_s + 0.25*self.peak_gain + 5*math.log10(self.z))

        # one-half the 3 dB beamwidth in the plane of interest
        self.psi_b = param.antenna_3_dB/2

        if self.l_s == -15:
            self.a = 2.58*math.sqrt(1 - 1.4*math.log10(self.z))
        elif self.l_s == -20:
            self.a = 2.58*math.sqrt(1 - 1.0*math.log10(self.z))
        elif self.l_s == -25:
            self.a = 2.58*math.sqrt(1 - 0.6*math.log10(self.z))
        elif self.l_s == -30:
            self.a = 2.58*math.sqrt(1 - 0.4*math.log10(self.z))
        else:
            sys.stderr.write("ERROR\nInvalid AntennaS1528 L_s parameter: " + str(self.l_s))
            sys.exit(1)

        self.b = 6.32
        self.alpha = 1.5

        self.x = self.peak_gain + self.l_s + 25*math.log10(self.b * self.psi_b)
        self.y = self.b * self.psi_b * math.pow(10, 0.04 * (self.peak_gain + self.l_s - self.l_f))

    def calculate_gain(self, *args, **kwargs) -> np.array:
        psi = np.absolute(kwargs["off_axis_angle_vec"])

        gain = np.zeros(len(psi))

        idx_0 = np.where(psi < self.a * self.psi_b)[0]
        gain[idx_0] = self.peak_gain - 3 * np.power(psi[idx_0] / self.psi_b, self.alpha)

        idx_1 = np.where((self.a * self.psi_b < psi) & (psi <= 0.5 * self.b * self.psi_b))[0]
        gain[idx_1] = self.peak_gain + self.l_s + 20 * math.log10(self.z)

        idx_2 = np.where((0.5 * self.b * self.psi_b < psi) & (psi <= self.b * self.psi_b))[0]
        gain[idx_2] = self.peak_gain + self.l_s

        idx_3 = np.where((self.b * self.psi_b < psi) & (psi <= self.y))[0]
        gain[idx_3] = self.x - 25 * np.log10(psi[idx_3])

        idx_4 = np.where((self.y < psi) & (psi <= 90))[0]
        gain[idx_4] = self.l_f

        idx_5 = np.where((90 < psi) & (psi <= 180))[0]
        gain[idx_5] = self.l_b

        return gain


if __name__ == '__main__':
    import matplotlib.pyplot as plt

    # initialize antenna parameters
    param = ParametersFssSs()
    param.antenna_gain = 39
    param.antenna_pattern = "ITU-R S.1528-0"
    param.antenna_3_dB = 2
    psi = np.linspace(0, 30, num = 1000)

    param.antenna_l_s = -15
    antenna = AntennaS1528(param)
    gain15 = antenna.calculate_gain(phi_vec=psi)

    param.antenna_l_s = -20
    antenna = AntennaS1528(param)
    gain20 = antenna.calculate_gain(phi_vec=psi)

    param.antenna_l_s = -25
    antenna = AntennaS1528(param)
    gain25 = antenna.calculate_gain(phi_vec=psi)

    param.antenna_l_s = -30
    antenna = AntennaS1528(param)
    gain30 = antenna.calculate_gain(phi_vec=psi)

    fig = plt.figure(figsize=(8,7), facecolor='w', edgecolor='k')  # create a figure object

    plt.plot(psi, gain15 - param.antenna_gain, "-b", label="$L_S = -15$ dB")
    plt.plot(psi, gain20 - param.antenna_gain, "-r", label="$L_S = -20$ dB")
    plt.plot(psi, gain25 - param.antenna_gain, "-g", label="$L_S = -25$ dB")
    plt.plot(psi, gain30 - param.antenna_gain, "-k", label="$L_S = -30$ dB")

    plt.ylim((-40, 10))
    plt.xlim((0, 30))
    plt.title("ITU-R S.1528-0 antenna radiation pattern")
    plt.xlabel("Relative off-axis angle, $\psi/\psi_{3dB}$")
    plt.ylabel("Gain relative to $G_{max}$ [dB]")
    plt.legend(loc="upper right")

#    ax = plt.gca()
#    ax.set_yticks([-30, -20, -10, 0])
#    ax.set_xticks(np.linspace(1, 9, 9).tolist() + np.linspace(10, 100, 10).tolist())

    plt.grid()
    plt.show()
