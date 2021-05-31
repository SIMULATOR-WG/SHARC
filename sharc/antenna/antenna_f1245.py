# -*- coding: utf-8 -*-
"""
@Created: Luciano Camilo Tue Fev 04 16:27:25 2021

"""

import numpy as np
from sharc.parameters.parameters_fs import ParametersFs
import math

class AntennaF1245(object):
    """
    Implements reference radiation patterns for fixed wireless system antennas for use in coordination studies use in
    interference assessment in the frequency range from 1 GHz to 86 GHz. (ITU-R F.1245-3)

    Attributes
    ----------
        antenna_gain (float): maximum gain of FS parabolic dish antenna
    """

    def __init__(self, param: ParametersFs):
        """
        Parameters
        ---------
            param (ParametersFS): antenna FS parameters
            phi : off-axis angle [deg]
        """
        super().__init__()
        self.peak_gain = param.antenna_gain
        lmbda = 3e8 / (param.frequency * 1e6)
        self.d_lmbda = param.diameter / lmbda

        self.g_l = 2 + 15 * math.log10(self.d_lmbda)
        self.phi_m = 20 / self.d_lmbda * math.sqrt(self.peak_gain - self.g_l)
        self.phi_r = 12.02 * math.pow(self.d_lmbda, -0.6)

    def calculate_gain(self, *args, **kwargs) -> np.array:
        phi = np.absolute(kwargs["off_axis_angle_vec"])

        if self.d_lmbda > 100:
            gain = self.calculate_gain_greater(phi)
        else:
            gain = self.calculate_gain_less(phi)

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
        gain[idx_0] = self.peak_gain - 0.0025 * np.power(self.d_lmbda * phi[idx_0], 2)

        idx_1 = np.where((self.phi_m <= phi) & (phi < self.phi_r))[0]
        gain[idx_1] = self.g_l

        idx_2 = np.where((self.phi_r <= phi) & (phi < 48))[0]
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
        gain[idx_1] = 39 - 5*math.log10(self.d_lmbda)-25 * np.log10(phi[idx_1])

        idx_2 = np.where((48 <= phi) & (phi <= 180))[0]
        gain[idx_2] = -3 - 5 * math.log10(self.d_lmbda)

        return gain


if __name__ == '__main__':
    import matplotlib.pyplot as plt

    phi = np.linspace(0.1, 180, num=100000)

    # initialize antenna parameters
    param_gt = ParametersFs()
    param_gt.antenna_pattern = "ITU-R F.699"
    param_gt.frequency = 1725
    param_gt.antenna_gain = 33
    param_gt.diameter = 3.2
    antenna_gt = AntennaF1245(param_gt)

    gain_gt = antenna_gt.calculate_gain(off_axis_angle_vec=phi)

    param_lt = ParametersFs()
    param_lt.antenna_pattern = "ITU-R F.1245"
    param_lt.frequency = 2600
    param_lt.antenna_gain = 25
    param_lt.diameter = 0.9
    antenna_lt = AntennaF1245(param_lt)
    gain_lt = antenna_lt.calculate_gain(off_axis_angle_vec=phi)
    csfont = {'fontname': 'Times New Roman'}
    hfont = {'fontname': 'Times New Roman'}
    # plt.legend(prop={'family': 'Times New Roman'})
    fig = plt.figure(facecolor='w', edgecolor='k')  # create a figure object

    #plt.semilogx(phi, gain_gt, "-b", label = "$f = 10.7$ $GHz,$ $D = 3$ $m$")
    plt.semilogx(phi, gain_lt, "--b", label="$f = 2032$ $MHz,$ $D = 4$ $m$")

    # plt.title("ITU-R F.699 antenna radiation pattern")
    plt.xlabel("Off-axis angle $\phi$ [degrees]", fontsize=18, color='black', **csfont)
    plt.ylabel("Gain relative to $G_m$ [dB]", fontsize=18, color='black', **csfont)
    # plt.legend(loc="lower left")
    plt.xlim((phi[0], phi[-1]))
    plt.ylim((-20, 40))

    pontos = 18000
    delta_theta = 180 / pontos
    theta = np.zeros(pontos)
    for i in range(pontos):
        theta[i] = i * delta_theta
    g = np.zeros(pontos)
    omega = np.zeros(pontos)
    for i in range(pontos):
        g[i] = 10 ** (gain_lt[i] / 10)

    omega = 2 * np.pi * np.sum(g * np.sin(theta * np.pi / 180)) * (delta_theta * np.pi / 180)
    print(omega)
    eff = (omega / (4 * np.pi))
    print(np.sum(eff))

    # ax = plt.gca()
    # ax.set_yticks([-30, -20, -10, 0])
    # ax.set_xticks(np.linspace(1, 9, 9).tolist() + np.linspace(10, 100, 10).tolist())

    plt.grid()
    plt.show()
