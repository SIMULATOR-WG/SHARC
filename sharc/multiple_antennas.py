# -*- coding: utf-8 -*-
"""
@Created: Luciano Camilo Tue Fev 04 16:27:25 2021

"""

import numpy as np
from sharc.parameters.parameters_fs import ParametersFs
from scipy import special
from sharc.antenna.antenna import Antenna
from sharc.parameters.parameters_fs import ParametersFs
import numpy as np
import math
import pandas as pd
import matplotlib.pyplot as plt

class AntennaF1245(object):
    """
    Implements a parabolic dish antenna pattern using Bessel of an FS antenna array following book Interference Analysis
    [John Paul] - Please there is an errata in this book on equation 3.79 and the correct equation is:
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


class AntennaF699(Antenna):
    """
    Implements reference radiation patterns for fixed wireless system antennas
    for use in coordination studies and interference assessment in the
    frequency range from 100 MHz to about 70 GHz. (ITU-R F.699-7)
    """

    def __init__(self, param: ParametersFs):
        super().__init__()
        self.peak_gain = param.antenna_gain
        lmbda = 3e8 / (param.frequency * 1e6)
        self.d_lmbda = param.diameter / lmbda

        self.g_l = 2 + 15 * math.log10(self.d_lmbda)
        self.phi_m = 20 / self.d_lmbda * math.sqrt(self.peak_gain - self.g_l)
        self.phi_r = 15.85 * math.pow(self.d_lmbda, -0.6)

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
        gain[idx_2] = 32 - 25 * np.log10(phi[idx_2])

        idx_3 = np.where((48 <= phi) & (phi <= 180))[0]
        gain[idx_3] = -10

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

        idx_1 = np.where((self.phi_m <= phi) & (phi < 100 / self.d_lmbda))[0]
        gain[idx_1] = self.g_l

        idx_2 = np.where((100 / self.d_lmbda <= phi) & (phi < 48))[0]
        gain[idx_2] = 52 - 10 * math.log10(self.d_lmbda) - 25 * np.log10(phi[idx_2])

        idx_3 = np.where((48 <= phi) & (phi <= 180))[0]
        gain[idx_3] = 10 - 10 * math.log10(self.d_lmbda)

        return gain


class AntennaBessel(object):
    """
    Implements a parabolic dish antenna pattern using Bessel of an FS antenna array following book Interference Analysis
    [John Paul] - Please there is an errata in this book on equation 3.79 and the correct equation is:
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
        #x = ((np.pi * self.d) / lambda1) * (np.sin(np.radians(phi)))
        x = ((np.pi * self.d) / lambda1) * (np.sin(np.radians(phi)))
        g = ((np.pi * self.d) / lambda1)**2
        g = 10 * np.log10(g)

        gain = 10 * np.log10((1 + np.cos(np.radians(phi))) * (special.j1(x) / x) ** 2)
        #gain = 10 * np.log10((2*(special.j1(x) / x) ** 2))
        #gain = 10 * np.log10(2*(special.j1(x) / x) ** 2)
        gain = gain + g
        return gain

        # phi = np.absolute(kwargs["off_axis_angle_vec"])
        # lambda1 = self.speed/(self.frequency*1E6)
        # x = ((np.pi*self.d)/lambda1)*(np.sin(np.radians(phi)))
        # gain = 10*np.log10((2*special.j1(x)/x)**2)
        # gain = gain + self.g_max
        # return gain


if __name__ == '__main__':

    from sharc.parameters.parameters_fs import ParametersFs
    from matplotlib import pyplot as plt

    data = pd.read_csv('medido2.txt', skiprows=0, sep='		', header=None)

    # plt.figure(figsize=(8, 6))

    x = data[0]
    y = data[1]
    print(x)
    phi = np.arange(-180, 180.5, step=0.5)
    print(len(y))
    print(len(phi))
    plt.grid(which='minor', alpha=0.5)
    plt.grid(which='major', alpha=0.5)
    plt.grid(True, color='b', linestyle='--', linewidth=0.2)
    plt.plot(phi, y, 'r-', linewidth=1.5, color='darkorange', label='Measured ')

    """
    Test routine - Comparison between SHARC code and Figure 3.36 of Interference Analysis book.
    """
    param = ParametersFs()
    param.antenna_gain = 35
    param.diameter = 4
    param.frequency = 2032
    #param.diameter = 2.4
    #param.frequency = 3600
    antenna = AntennaBessel(param)
    phi = np.arange(-180, 180, step=0.1)
    #phi = np.arange(0, 90, step=0.1)
    csfont = {'fontname': 'Times New Roman'}
    hfont = {'fontname': 'Times New Roman'}
    plt.figure(1)
    # plt.legend(loc='upper left')
    plt.plot(phi, antenna.calculate_gain(off_axis_angle_vec=phi), '--', color='black',linewidth=0.4,
             label='Bessel')

    #plt.semilogx(phi, antenna.calculate_gain(off_axis_angle_vec=phi), 'b--',
    #             label="Bessel $f = 2032$ $MHz,$ $D = 4$ $m$")
    plt.legend(loc="lower left")
    plt.xlim(-180, 180)
    #plt.xlim(0, 90)
    plt.ylim(-50, 40)
    # plt.xlim(-90, 90)
    plt.grid(which='minor', alpha=0.2)
    plt.grid(which='major', alpha=0.5)
    plt.grid(True, color='black', linestyle='--', linewidth=0.1)
    plt.title('Comparison Between Antenna Patterns')
    plt.xlabel("Off-axis angle $\phi$ [degrees]", fontsize=16, color='black', **csfont)
    plt.ylabel("Gain relative to $G_m$ [dB]", fontsize=16, color='black', **csfont)
    # plt.legend()


   # phi = np.linspace(-180, 180, num=100000)



    param_lt = ParametersFs()
    param_lt.antenna_pattern = "ITU-R F.699"
    param_lt.frequency = 2032
    param_lt.antenna_gain = 35.6
    param_lt.diameter = 4
    antenna_lt = AntennaF699(param_lt)
    gain_lt = antenna_lt.calculate_gain(off_axis_angle_vec=phi)
    csfont = {'fontname': 'Times New Roman'}
    hfont = {'fontname': 'Times New Roman'}
    #plt.legend(prop={'family': 'Times New Roman'})
    #fig = plt.figure( facecolor='w', edgecolor='k')  # create a figure object

    #plt.semilogx(phi, gain_gt, "--b", label = "$f = 10.7$ $GHz,$ $D = 3$ $m$")
    #plt.semilogx(phi, gain_lt, "--r", label="ITU-R F.699 $f = 2032$ $MHz,$ $D = 4$ $m$")
    #plt.plot(phi, gain_lt, "--r", label="ITU-R F.699 $f = 2032$ $MHz,$ $D = 4$ $m$")
    plt.plot(phi, gain_lt, "--",color='royalblue', linewidth=1.5, label="ITU-R F.699")

    # plt.title("ITU-R F.699 antenna radiation pattern")
    plt.xlabel("Off-axis angle $\phi$ [degrees]",fontsize=16, color='black', **csfont)
    plt.ylabel("Gain relative to $G_m$ [dB]",fontsize=16, color='black', **csfont)
    plt.legend(loc="lower left")
    #plt.xlim((phi[0], phi[-1]))
    #plt.ylim((-10, 40))
    plt.xlim(-180, 180)
    #plt.xlim(0, 90)
    plt.ylim(-50, 40)
    # ax = plt.gca()
    # ax.set_yticks([-30, -20, -10, 0])
    # ax.set_xticks(np.linspace(1, 9, 9).tolist() + np.linspace(10, 100, 10).tolist())

    antenna_1245 = AntennaF1245(param_lt)
    gain_f1245 = antenna_1245.calculate_gain(off_axis_angle_vec=phi)
    csfont = {'fontname': 'Times New Roman'}
    hfont = {'fontname': 'Times New Roman'}
    # plt.legend(prop={'family': 'Times New Roman'})
    #fig = plt.figure(facecolor='w', edgecolor='k')  # create a figure object

    # plt.semilogx(phi, gain_gt, "-b", label = "$f = 10.7$ $GHz,$ $D = 3$ $m$")
    plt.plot(phi, gain_f1245, "--b", color='red', label="ITU-R F.1245")
    plt.legend(loc="lower left")
    plt.grid()


    plt.show()
