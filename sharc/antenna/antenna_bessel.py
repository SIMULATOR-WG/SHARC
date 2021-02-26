# -*- coding: utf-8 -*-
"""
@Created: Luciano Camilo Tue Fev 04 16:27:25 2021

"""

import numpy as np
from sharc.parameters.parameters_fs import ParametersFs
from scipy import special
import math

class AntennaBessel(object):
    """
    Implements a parabolic dish antenna pattern, Bessel model, based on Interference Analysis[John Pahl] book
    - Please there is an errata in this book on equation 3.79 and the correct equation is:
    G(theta) = 10*log10((2j1(x)/x)^2)  - Equation 3.79

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
        self.g_max = param.antenna_gain
        self.d = param.diameter
        self.frequency = param.frequency
        self.speed = 299792458

    def calculate_gain(self, **kwargs) -> np.array:

        phi = np.absolute(kwargs["off_axis_angle_vec"])
        lambda1 = self.speed / (self.frequency * 1E6)
        x = ((np.pi * self.d) / lambda1) * (np.sin(np.radians(phi)))
        g = ((np.pi * self.d) / lambda1) ** 2
        g = 10 * np.log10(g)
        gain = 10 * np.log10((1 + np.cos(np.radians(phi))) * (special.j1(x) / x) ** 2)
        gain = gain + g
        return gain


if __name__ == '__main__':

    from sharc.parameters.parameters_fs import ParametersFs
    from matplotlib import pyplot as plt
    """
    Test routine - Comparison between SHARC code and Figure 3.36 of Interference Analysis book.
    """

    param = ParametersFs()
    param.antenna_gain = 35
    param.diameter = 4
    param.frequency = 2032
    antenna = AntennaBessel(param)
    phi = np.arange(0, 180, step=0.05)
    csfont = {'fontname': 'Times New Roman'}
    hfont = {'fontname': 'Times New Roman'}
    plt.figure(1)
    # plt.legend(loc='upper left')
    #plt.plot(phi, antenna.calculate_gain(off_axis_angle_vec=phi), 'b--',
    #         label='Bessel')

    plt.semilogx(phi, antenna.calculate_gain(off_axis_angle_vec=phi), 'b--',
                 label="Bessel $f = 2032$ $MHz,$ $D = 4$ $m$")
    plt.legend(loc="lower left")
    plt.xlim(0, 180)
    plt.ylim(-50, 40)
    # plt.xlim(-90, 90)
    plt.grid(which='minor', alpha=0.2)
    plt.grid(which='major', alpha=0.5)
    plt.grid(True, color='b', linestyle='--', linewidth=0.2)
    plt.title('Theoretical Bessel Antenna Pattern')
    plt.xlabel("Off-axis angle $\phi$ [degrees]", fontsize=16, color='black', **csfont)
    plt.ylabel("Gain relative to $G_m$ [dB]", fontsize=16, color='black', **csfont)
    # plt.legend()
    plt.show()
