# -*- coding: utf-8 -*-
"""
@Created: Luciano Camilo on Wed Marc 24 19:58:58 2021

"""

import numpy as np
from sharc.parameters.parameters_arns import ParametersArns


class AntennaCosineRadar(object):
    """
    Implements a cosine antenna pattern of an Aeronautical Surveillance Radar antenna following
    ITU-R M.1851-1, Table 2/3 and Figure 8.

    Attributes
    ----------
        antenna_gain (float): maximum gain of Radar antenna
    """

    def __init__(self, param: ParametersArns):
        """
        Constructs an AntennaElement object.

        Parameters
        ---------
            param (ParametersArns): antenna ARNS parameters
        """
        super().__init__()
        self.g_max = param.antenna_gain
        self.phi_3db = param.beamwidth
        self.maximum_csc = param.csc2_angle

    def calculate_gain(self, **kwargs) -> np.array:

        phi = np.absolute(kwargs["off_axis_angle_vec"])
        const = np.pi * 68.8
        g = np.zeros(len(phi))
        mi = np.zeros(len(phi))

        for i in range(len(phi)):

            if -self.phi_3db <= phi[i] <= self.phi_3db:
                mi[i] = const*np.sin(np.radians(phi[i]))/self.phi_3db
                g[i] = 20*np.log10((np.pi/2)*((np.cos(mi[i]))/((np.pi/2)**2-(mi[i])**2)))+4.32-0.392
            if self.phi_3db <= phi[i] <= 180 or -self.phi_3db >= phi[i] >= -180:
                g[i] = -17.51*np.log(2.33*(np.abs((phi[i]))/self.phi_3db))-4.32  # cos2
                if g[i] < -50:
                    g[i] = -50
        g = g+self.g_max

        return g


if __name__ == '__main__':

    from sharc.parameters.parameters_arns import ParametersArns
    from matplotlib import pyplot as plt
    """
    Test routine - Comparison between SHARC code and Figure 8 of ITU R - M.1851-1
    """
    param = ParametersArns()
    param.antenna_gain = 0
    param.beamwidth = 3
    param.csc2_angle = 30
    antenna = AntennaCosineRadar(param)
    phi = np.arange(-180, 180, step=0.001)
    # theta = np.arange(-90, 90, step=0.1)
    csfont = {'fontname': 'Times New Roman'}
    hfont = {'fontname': 'Times New Roman'}
    plt.figure(1)

    # plt.legend(loc='upper left')
    # plt.plot(theta, antenna.calculate_gain(theta_vec=theta), 'b--', color = 'darkorange',
    # label='$CSC^2$ Antenna Pattern')
    plt.plot(phi, antenna.calculate_gain(off_axis_angle_vec=phi), 'b-', color='darkorange',
             label='$COS$ Antenna Pattern')
    # plt.semilogx(phi, antenna.calculate_gain(off_axis_angle_vec=phi), 'b--', color = 'darkorange',
    #  label='$CSC^2$ Antenna Pattern')
    # plt.xlim(-180, 180)
    # plt.ylim(-30, 0)
    plt.xlim(0, 30)
    plt.grid(which='minor', alpha=0.5)
    plt.grid(which='major', alpha=0.5)
    plt.grid(True, color='b', linestyle='--', linewidth=0.2)
    plt.title('Aeronautical Surveillance Radar COS Antenna Pattern')
    plt.xlabel('Azimuth (degrees)', fontsize=16, color='black', **csfont)
    plt.ylabel('Gain (dBi)', fontsize=16, color='black', **csfont)
    plt.legend()
    plt.show()
