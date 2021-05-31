# -*- coding: utf-8 -*-
"""
@Created: Luciano Camilo on Fri May 28 15:54:00 2021
"""
from sharc.antenna.antenna import Antenna
from sharc.parameters.parameters_arns import ParametersArns

import numpy as np


class AntennaRadarPhasedArray(object):
    """
    Implements a Phased array antenna pattern of an Aeronautical Surveillance Radar antenna following
    ITU-R M.1851-1, Equation 20 and Figure 21.

    Attributes
    ----------
        antenna_gain (float): maximum gain of Radar antenna
    """

    def __init__(self, param: ParametersArns):
        """
        Constructs an AntennaRadar object.
        Does not receive angles in local coordinate system.
        Elevation taken wrt x-y plane.

        Parameters
        ---------
            azimuth (float): antenna's physical azimuth inclination
            elevation (float): antenna's physical elevation inclination
        """
        super().__init__()

        self.g_max = param.antenna_gain
        if self.g_max == 41:
            self.normalization = 0.833
            self.number_elements = 38
        if self.g_max == 46:
            self.normalization = 3.48681
            self.number_elements = 70
        self.theta_3db = 1.1
        self.phi_3db = 1.1
        self.frequency = param.frequency
        self.l_ambda = 299792458 / (self.frequency * 1e6)
        self.element_space = param.element_space * self.l_ambda
        self.beamsteeringangle_az = param.beamsteeringangle_az
        self.beamsteeringangle_el = param.beamsteeringangle_el

    def calculate_gain(self, *args, **kwargs) -> np.array:
        """
        Calculates the gain in the given direction.

        Parameters
        ----------
            phi_vec (np.array): azimuth angles [degrees]
            theta_vec (np.array): elevation angles [degrees]

        Returns
        -------
            gain (np.array): gain corresponding to each of the given directions
        """
        phi = np.asarray(kwargs["off_axis_angle_vec"])
        theta = np.asarray(kwargs["theta_vec"])

        phi_l, theta_l = self.to_local_coord(phi, theta)
        gain = self.g_max + self.get_array_gain_az(phi_l) + self.get_array_gain_el(phi_l)

        # print(theta_l)
        # print(self.g_max + self.get_gain_elev(phi_l))
        return gain

    def get_array_gain_az(self, theta: np.array) -> np.array:
        """
        Returns the antenna gain in the azimuth plane for the given direction.
        """
        const = np.pi * 50.8
        g = np.zeros(len(theta))
        mi = np.zeros(len(theta))
        af = np.zeros(len(theta))
        psi = np.zeros(len(theta))

        for i in range(len(theta)):
            mi[i] = (const * np.sin(np.radians(theta[i]))) / self.theta_3db
            g[i] = 20 * np.log10(((np.sin(np.radians(mi[i]))) / (mi[i])))  # uniform
            psi[i] = (2 * np.pi * (self.element_space / self.l_ambda) * (
                    np.sin((np.radians(theta[i]))) - np.sin(np.radians(self.beamsteeringangle_az))))
            af[i] = np.sin(((self.number_elements * psi[i]) / 2)) / np.sin((psi[i] / 2))
            g[i] = g[i] + 10 * np.log10((1 / self.number_elements)) + 10 * np.log10(abs(af[i]) ** 2) + 20.2 - self.normalization

        return g

    def get_array_gain_el(self, theta: np.array) -> np.array:
        """
        Returns the antenna gain in the elevation plane for the given direction.
        """

        const = np.pi * 50.8
        g = np.zeros(len(theta))
        mi = np.zeros(len(theta))
        af = np.zeros(len(theta))
        psi = np.zeros(len(theta))

        for i in range(len(theta)):
            mi[i] = (const * np.sin(np.radians(theta[i]))) / self.theta_3db
            g[i] = 20 * np.log10(((np.sin(np.radians(mi[i]))) / (mi[i])))  # uniform
            psi[i] = (2 * np.pi * (self.element_space / self.l_ambda) * (
                np.sin((np.radians(theta[i]))) - np.sin(np.radians(self.beamsteeringangle_el))))
            af[i] = np.sin(((self.number_elements * psi[i]) / 2)) / np.sin((psi[i] / 2))
            g[i] = g[i] + 10 * np.log10((1 / self.number_elements)) + 10 * np.log10(abs(af[i]) ** 2) + 20.2 - self.normalization

        return g

    def to_local_coord(self, phi: np.array, theta: np.array) -> tuple:
        """
        Converts the azimuth and elevation angles to the local coordinate system,
        taking into account the physical orientation of the antenna.
        """

        # phi_l = phi
        # theta_l = theta-90
        # theta_l = -theta

        phi_l = phi
        theta_l = -theta

        return phi_l, theta_l


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation

    param = ParametersArns()
    param.antenna_gain = 41
    param.beamwidth_el = 1.1
    param.beamwidth_az = 1.1
    param.frequency = 2700
    param.beamsteeringangle_az = 0
    param.beamsteeringangle_el = 0
    #param.number_elements = 70
    param.number_elements = 38
    param.element_space = 0.5

    azimuth = np.linspace(-90, 90, num=5000)
    elevation = np.linspace(-90, 90, num=5000)

    gain = 0
    azimuth_l = 0
    elevation_l = 0

    antenna = AntennaRadarPhasedArray(param)

    gain_az = gain + antenna.get_array_gain_az(azimuth)
    gain_elev = gain + antenna.get_array_gain_el(elevation)

    fig = plt.figure(figsize=(15, 5), facecolor='w', edgecolor='k')

    ax1 = fig.add_subplot(121)
    #ax1.plot(azimuth, gain_az, 'b--', color='blue', label='Phased Array Radar Antenna')
    ax1.semilogx(azimuth, gain_az, 'b--', color='blue', label='Phased Array Radar Antenna')
    plt.legend(loc='upper right')
    ax1.grid(True)
    ax1.set_xlabel("Azimuth angle [deg]")
    ax1.set_ylabel("Antenna gain [dBi]")
    plt.title("Phased Array Radar Antenna")
    ax1.set_xlim([-90, 90])
    #ax1.set_xlim([0, 90])
    ax1.set_ylim([-50, 10])

    ax2 = fig.add_subplot(122)
    #ax2.plot(elevation, gain_elev, 'b--', color='darkorange', label='Phased Array Radar Antenna')
    ax2.semilogx(elevation, gain_elev, 'b--', color='darkorange', label='Phased Array Radar Antenna')
    ax2.grid(True)
    ax2.set_xlabel("Elevation angle [deg]")
    ax2.set_ylabel("Antenna gain [dBi]")
    ax2.set_xlim([0, 90])
    ax2.set_ylim([-50, 10])
    plt.title("Phased Array Radar Antenna")
    plt.legend(loc='upper right')
    plt.show()
    # i=0
    # plt.ion()
    # import time
    # while i<80:
    #     param.beamsteeringangle_az = i
    #     antenna = AntennaRadarPhasedArray(param)
    #     gain_az = gain + antenna.get_array_gain_az(azimuth)
    #     plt.plot(azimuth, gain_az, 'b--', color='blue', label='Phased Array Radar Antenna')
    #     #plt.show()
    #     #plt.draw()
    #     plt.draw()
    #     plt.pause(0.001)
    #     plt.clf()
    #     plt.xlim([-90, 90])
    #     plt.ylim([-60, 0])
    #     #ax1.clear()
    #     #plt.pause(0.1)
    #     #plt.clf()
    #     #plt.clf()
    #     #plt.close()
    #     i += 1

# plt.show()
