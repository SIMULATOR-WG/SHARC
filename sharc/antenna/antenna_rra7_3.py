# -*- coding: utf-8 -*-

from sharc.antenna.antenna import Antenna
from sharc.parameters.parameters_antenna_with_diameter import ParametersAntennaWithDiameter
from sharc.antenna.antenna_s465 import AntennaS465

import numpy as np
import math


class AntennaReg_RR_A7_3(Antenna):
    """
    Implements the Earth station antenna pattern in the MetSat service
    according to Recommendation ITU Radio Regulations Appendix 7, Annex 3
    """

    def __init__(self, param: ParametersAntennaWithDiameter):
        super().__init__()
        self.peak_gain = param.antenna_gain
        lmbda = 3e8 / (param.frequency * 1e6)
        if param.diameter:
            D_lmbda = param.diameter / lmbda
        else:
            # From Note 1, we can estimate D_lmbda from
            # 20 log(D_lmbda) =(aprox.)= G_max - 7.7
            D_lmbda = math.pow(10, ((self.peak_gain - 7.7) / 20))

        self.D_lmbda = D_lmbda

        if D_lmbda >= 100:
            self.g1 = -1 + 15 * np.log10(D_lmbda)
            self.phi_r = 15.85 * math.pow(D_lmbda, -0.6)
        elif D_lmbda >= 35:
            self.g1 = -21 + 25 * np.log10(D_lmbda)
            self.phi_r = 100 / D_lmbda
        else:
            raise ValueError(f"Recommendation does not define antenna pattern when D/lmbda = {D_lmbda}")

        self.phi_m = 20 / D_lmbda * np.sqrt(self.peak_gain - self.g1)

        # if np.sqrt(self.peak_gain - self.g1) >= 5, then phi_m >= phi_r, and that may be a problem
        # if this is erroring in your simulation, you should check with a professor how to deal with this
        # since the document doesn't specify
        if self.phi_m >= self.phi_r:
            raise ValueError(f"Recommendation doesn't specify what to do when phi_m ({self.phi_m}) >= phi_r ({self.phi_r})")

    def calculate_gain(self, *args, **kwargs) -> np.array:
        phi = np.absolute(kwargs["off_axis_angle_vec"])

        gain = np.zeros(phi.shape)

        idx_00 = np.where(phi < self.phi_m)[0]
        gain[idx_00] = self.peak_gain - 2.5 * self.D_lmbda * self.D_lmbda * phi[idx_00] * phi[idx_00] / 1000

        idx_0 = np.where((self.phi_m <= phi) & (phi < self.phi_r))[0]
        gain[idx_0] = self.g1

        idx_1 = np.where((self.phi_r <= phi) & (phi < 36))[0]
        gain[idx_1] = 29 - 25 * np.log10(phi[idx_1])

        idx_2 = np.where((36 <= phi) & (phi <= 180))[0]
        gain[idx_2] = -10

        return gain


if __name__ == '__main__':
    import matplotlib.pyplot as plt

    phi = np.linspace(0.1, 100, num=100000)

    # initialize antenna parameters
    param27 = ParametersAntennaWithDiameter()
    param27.antenna_pattern = "ITU-R S.465-6"
    param27.frequency = 7500
    param27.antenna_gain = 50
    param27.diameter = 5
    antenna27 = AntennaS465(param27)

    gain27 = antenna27.calculate_gain(off_axis_angle_vec=phi)

    param43 = ParametersAntennaWithDiameter()
    param43.antenna_pattern = "ITU-R S.465-6"
    param43.frequency = param27.frequency
    param43.antenna_gain = 59
    param43.diameter = 13
    antenna43 = AntennaReg_RR_A7_3(param43)
    gain43 = antenna43.calculate_gain(off_axis_angle_vec=phi)

    fig = plt.figure(
        figsize=(8, 7), facecolor='w',
        edgecolor='k',
    )  # create a figure object

    plt.semilogx(
        phi, gain27, "-b",
        label="ES for Sat Q (ITU-R S.465-6)",
    )
    plt.semilogx(
        phi, gain43, "-r",
        label="ES for Sat P (ITU-R Reg. R.R. Appendix 7, Annex 3)",
    )

    plt.title("ES antenna radiation patterns")
    plt.xlabel(r"Off-axis angle $\phi$ [deg]")
    plt.ylabel("Gain [dBi]")
    plt.legend(loc="lower left")
    plt.xlim((phi[0], phi[-1]))
    # plt.ylim((-80, 10))

    # ax = plt.gca()
    # ax.set_yticks([-30, -20, -10, 0])
    # ax.set_xticks(np.linspace(1, 9, 9).tolist() + np.linspace(10, 100, 10).tolist())

    plt.grid()
    plt.show()
