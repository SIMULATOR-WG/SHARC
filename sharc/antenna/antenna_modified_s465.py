# -*- coding: utf-8 -*-
"""
Created on Thu Aug 24 18:41:17 2017

@author: edgar
"""

from sharc.antenna.antenna import Antenna
from sharc.parameters.parameters_fss_es import ParametersFssEs

import numpy as np
import math

class AntennaModifiedS465(Antenna):
    """
    Implements the Earth station antenna pattern in the fixed-satellite service
    according to Recommendation ITU-R S.465-6 Annex 1
    """

    def __init__(self, param: ParametersFssEs):
        super().__init__()
        self.peak_gain = param.antenna_gain
        self.phi_min = 5.0 # degrees
        self.envelope_gain = param.antenna_envelope_gain
        self.envelope_angle = 10**((32 - self.envelope_gain) / 25.)

    def calculate_gain(self, *args, **kwargs) -> np.array:

        phi = np.absolute(kwargs["off_axis_angle_vec"])

        gain = np.zeros(phi.shape)

        idx_0 = np.where(phi < self.phi_min)[0]
        gain[idx_0] = self.peak_gain

        idx_1 = np.where((self.phi_min <= phi) & (phi < self.envelope_angle))[0]
        gain[idx_1] = 32 - 25 * np.log10(phi[idx_1])

        idx_2 = np.where((self.envelope_angle <= phi) & (phi <= 180))[0]
        gain[idx_2] = self.envelope_gain

        return gain

if __name__ == '__main__':
    import matplotlib.pyplot as plt

    phi = np.linspace(0.1, 180, num = 100000)

    # initialize antenna parameters
    param0 = ParametersFssEs()
    param0.frequency = 3628
    param0.antenna_gain = 32
    param0.antenna_envelope_gain = 0
    antenna0 = AntennaModifiedS465(param0)
    gain0 = antenna0.calculate_gain(off_axis_angle_vec = phi)

    param4 = ParametersFssEs()
    param4.frequency = 3628
    param4.antenna_gain = 32
    param4.antenna_envelope_gain = -4
    antenna4 = AntennaModifiedS465(param4)
    gain4 = antenna4.calculate_gain(off_axis_angle_vec = phi)

    param10 = ParametersFssEs()
    param10.frequency = 3628
    param10.antenna_gain = 32
    param10.antenna_envelope_gain = -10
    antenna10 = AntennaModifiedS465(param10)
    gain10 = antenna10.calculate_gain(off_axis_angle_vec = phi)
    
    
    fig = plt.figure(figsize=(8,7), facecolor='w', edgecolor='k')  # create a figure object

    plt.semilogx(phi, gain0, "-b", label = "$G_{env}$ = 0 dBi")
    plt.semilogx(phi, gain4, "-r", label = "$G_{env}$ = -4 dBi")
    plt.semilogx(phi, gain10, "-g", label = "$G_{env}$ = -10 dBi")

    plt.title("Modified ITU-R S.465-6 antenna radiation pattern")
    plt.xlabel("Off-axis angle $\phi$ [deg]")
    plt.ylabel("Gain relative to $G_m$ [dB]")
    plt.legend(loc="upper right")
    plt.xlim((1, 180))
    plt.ylim((-20, 40))

    #ax = plt.gca()
    #ax.set_yticks([-30, -20, -10, 0])
    #ax.set_xticks(np.linspace(1, 9, 9).tolist() + np.linspace(10, 100, 10).tolist())

    plt.grid()
    plt.show()
