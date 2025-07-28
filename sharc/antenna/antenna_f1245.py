
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 4 17:08:00 2018
@author: Calil
"""
import matplotlib.pyplot as plt
from sharc.antenna.antenna import Antenna
from sharc.parameters.imt.parameters_imt import ParametersImt
from sharc.parameters.imt.parameters_antenna_imt import ParametersAntennaImt
import numpy as np
import math


class AntennaF1245(Antenna):
    """
    Implements reference radiation patterns for HAPS gateway antennas
    for use in coordination studies and interference assessment in the
    frequency range from 1 GHz to about 70 GHz. (ITU-R F.1245-2)
    """

    def __init__(self, param: ParametersImt, param_ant: ParametersAntennaImt):
        super().__init__()
        self.peak_gain = param_ant.peak_gain
        lmbda = 3e8 / (param.frequency * 1e6)
        self.d_lmbda = param_ant.diameter / lmbda
        self.g_l = 2 + 15 * math.log10(self.d_lmbda)
        self.phi_m = (20 / self.d_lmbda) * math.sqrt(self.peak_gain - self.g_l)
        self.phi_r = 12.02 * math.pow(self.d_lmbda, -0.6)

    def calculate_gain(self, *args, **kwargs) -> np.array:
        phi_vec = np.absolute(kwargs["phi_vec"])
        theta_vec = np.absolute(kwargs["theta_vec"])
        beams_l = np.absolute(kwargs["beams_l"])
        off_axis = self.calculate_off_axis_angle(phi_vec, theta_vec)
        if self.d_lmbda > 100:
            gain = self.calculate_gain_greater(off_axis)
        else:
            gain = self.calculate_gain_less(off_axis)
            idx_max_gain = np.where(beams_l == -1)[0]
            gain[idx_max_gain] = self.peak_gain
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
        gain[idx_0] = self.peak_gain - 2.5e-3 * \
            np.power(self.d_lmbda * phi[idx_0], 2)
        phi_thresh = max(self.phi_m, self.phi_r)
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
        gain[idx_0] = self.peak_gain - 0.0025 * \
            np.power(self.d_lmbda * phi[idx_0], 2)
        idx_1 = np.where((self.phi_m <= phi) & (phi < 48))[0]
        gain[idx_1] = 39 - 5 * \
            math.log10(self.d_lmbda) - 25 * np.log10(phi[idx_1])
        idx_2 = np.where((48 <= phi) & (phi < 180))[0]
        gain[idx_2] = -3 - 5 * math.log10(self.d_lmbda)
        return gain

    def add_beam(self, phi: float, theta: float):
        self.beams_list.append((phi, theta))

    def calculate_off_axis_angle(self, Az, b):
        Az0 = self.beams_list[0][0]
        a = 90 - self.beams_list[0][1]
        C = Az0 - Az
        off_axis_rad = np.arccos(
            np.cos(np.radians(a)) * np.cos(np.radians(b)) + np.sin(np.radians(a)) *
            np.sin(np.radians(b)) * np.cos(np.radians(C)),
        )
        off_axis_deg = np.degrees(off_axis_rad)
        return off_axis_deg


if __name__ == '__main__':
    phi = np.linspace(0.1, 180, num=100000)
    theta = 90 * np.ones_like(phi)
    beams_idx = np.zeros_like(phi, dtype=int)
    # initialize antenna parameters
    param = ParametersImt()
    param.frequency = 10700
    param_gt = ParametersAntennaImt()
    param_gt.peak_gain = 49.8
    param_gt.diameter = 3
    antenna_gt = AntennaF1245(param, param_gt)
    antenna_gt.add_beam(0, 0)
    gain_gt = antenna_gt.calculate_gain(
        phi_vec=phi,
        theta_vec=theta,
        beams_l=beams_idx,
    )
    param.frequency = 27500
    param_lt = ParametersAntennaImt()
    param_lt.peak_gain = 36.9
    param_lt.diameter = 0.3
    antenna_lt = AntennaF1245(param, param_lt)
    antenna_lt.add_beam(0, 0)
    gain_lt = antenna_lt.calculate_gain(
        phi_vec=phi,
        theta_vec=theta,
        beams_l=beams_idx,
    )
    fig = plt.figure(
        figsize=(8, 7), facecolor='w',
        edgecolor='k',
    )  # create a figure object
    plt.semilogx(phi, gain_gt, "-b", label="$f = 10.7$ $GHz,$ $D = 3$ $m$")
    plt.semilogx(phi, gain_lt, "-r", label="$f = 27.5$ $GHz,$ $D = 0.3$ $m$")
    plt.title("ITU-R F.1245 antenna radiation pattern")
    plt.xlabel(r"Off-axis angle $\phi$ [deg]")
    plt.ylabel("Gain relative to $G_m$ [dB]")
    plt.legend(loc="lower left")
    plt.xlim((phi[0], phi[-1]))
    plt.ylim((-20, 50))
    # ax = plt.gca()
    # ax.set_yticks([-30, -20, -10, 0])
    # ax.set_xticks(np.linspace(1, 9, 9).tolist() + np.linspace(10, 100, 10).tolist())
    plt.grid()
    plt.show()
