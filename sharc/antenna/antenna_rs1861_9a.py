# -*- coding: utf-8 -*-
"""
Created on Thu May  9 17:24:39 2019

@author: edgar
"""

from sharc.antenna.antenna import Antenna
from sharc.parameters.parameters_eess_passive import ParametersEessPassive

import numpy as np

class AntennaRS1861_9A(Antenna):
    """
    Implements the reference antenna pattern described in Figure 9a from
    Recommendation ITU-R RS.1861.     
    """

    def __init__(self, param: ParametersEessPassive):
        super().__init__()
        self.peak_gain = param.antenna_gain


    def calculate_gain(self, *args, **kwargs) -> np.array:
        phi = np.absolute(kwargs["off_axis_angle_vec"])

        gain = np.zeros(phi.shape)

        id0 = np.where(phi <= 2.5)[0]
        gain[id0] = self.peak_gain - 5.2*np.power(phi[id0], 2) - 0.6*(phi[id0])
        
        id1 = np.where((2.5 < phi) & (phi <= 10))[0]
        gain[id1] = self.peak_gain - 3.47*phi[id1] - 25.3
        
        id2 = np.where((phi > 10))[0]
        gain[id2] = self.peak_gain - 60
        
        return gain


if __name__ == '__main__':
    import matplotlib.pyplot as plt

    phi = np.linspace(0, 18, num = 10000)

    # initialize antenna parameters
    param = ParametersEessPassive()
    param.antenna_pattern = "ITU-R RS.1861 Fig 9a"
    param.antenna_gain = 40
    antenna = AntennaRS1861_9A(param)

    gain = antenna.calculate_gain(off_axis_angle_vec = phi)

    fig = plt.figure(figsize=(8,5), facecolor='w', edgecolor='k')  # create a figure object
    plt.plot(phi, gain - param.antenna_gain, "-b")

    plt.title("ITU-R RS.1861 Fig 9a antenna radiation pattern")
    plt.xlabel("Off-axis angle $\phi$ [deg]")
    plt.ylabel("Normalized antenna gain [dBi]")
    plt.xlim((phi[0], phi[-1]))
    plt.ylim((-60, 0))

    plt.grid()
    plt.show()
