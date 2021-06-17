# -*- coding: utf-8 -*-
"""
@Created: Luciano Camilo on Mon Apr 05 14:08:00 2021

"""

from sharc.antenna.antenna import Antenna
from sharc.parameters.parameters_arns import ParametersArns

import numpy as np


class AntennaMeteorologicalRadarUniform(Antenna):
    """
    Implements a Uniform Meteorological Radar antenna pattern following ITU-R M.1851-1, Table 5.
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
        self.theta_3db = param.beamwidth_el
        self.phi_3db = param.beamwidth_az

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
        gain = self.g_max + self.get_gain_az(phi_l)+self.get_gain_elev(phi_l)
        # print(theta_l)
        # print(self.g_max + self.get_gain_elev(phi_l))
        # print(gain)
        return gain

    def get_gain_az(self, phi: np.array) -> np.array:
        """
        Returns the antenna gain in the azimuth plane for the given direction.
        """
        const = np.pi * 50.8
        g = np.zeros(len(phi))
        mi = np.zeros(len(phi))

        for i in range(len(phi)):
            if -self.phi_3db <= phi[i] <= self.phi_3db-0.5:
                mi[i] = const * np.sin(np.radians(phi[i])) / np.radians(self.phi_3db)
                g[i] = 10 * np.log10((np.sin(np.radians(mi[i])) / mi[i]) ** 2) + 35.2
            if self.phi_3db-0.5 <= phi[i] <= 180 or -self.phi_3db+0.5 >= phi[i] >= -180:
                g[i] = -8.584 * np.log(2.876 * (np.abs((phi[i])) / self.phi_3db)) - 3.72  # cos2
                if g[i] < -30:
                    g[i] = -30
        return g

    def get_gain_elev(self, theta: np.array) -> np.array:
        """
        Returns the antenna gain in the elevation plane for the given direction.
        """

        const = np.pi * 50.8
        g = np.zeros(len(theta))
        mi = np.zeros(len(theta))

        for i in range(len(theta)):
            if -self.theta_3db <= theta[i] <= self.theta_3db-0.5:
                mi[i] = const * np.sin(np.radians(theta[i])) / np.radians(self.theta_3db)
                g[i] = 10 * np.log10((np.sin(np.radians(mi[i])) / mi[i]) ** 2) + 35.2
            if self.theta_3db-0.5 <= theta[i] <= 180 or -self.theta_3db+0.5 >= theta[i] >= -180:
                g[i] = -8.584 * np.log(2.876 * (np.abs((theta[i])) / self.theta_3db)) - 3.72  # cos2
                if g[i] < -30:
                    g[i] = -30
        return g

    def theoretical(self, theta: np.array) -> np.array:
        """
        Returns the theoretical antenna gain according Rec. M. 1851-2.
        """

        const = np.pi * 50.8
        g = np.zeros(len(theta))
        mi = np.zeros(len(theta))

        for i in range(len(theta)):
            mi[i] = const * np.sin(np.radians(theta[i])) / np.radians(self.theta_3db)
            g[i] = 10*np.log10((np.sin(np.radians(mi[i])) / mi[i])**2)+35.2
            # print(g[i])
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

    param = ParametersArns()
    param.antenna_gain = 0
    param.beamwidth_el = 3
    param.beamwidth_az = 3

    azimuth = np.linspace(-90, 90, num=15000)
    elevation = np.linspace(-90, 90, num=15000)

    gain = 0
    azimuth_l = 2
    elevation_l = 2

    antenna = AntennaMeteorologicalRadarUniform(param)

    gain_az = gain + antenna.get_gain_az(azimuth)
    gain_elev = gain + antenna.get_gain_elev(elevation)
    theoretical = gain + antenna.theoretical(elevation)

    fig = plt.figure(figsize=(15, 5), facecolor='w', edgecolor='k')

    ax1 = fig.add_subplot(121)
    ax1.plot(azimuth, gain_az, 'b--', color='blue', label='Uniform')
    ax1.plot(azimuth, theoretical, 'b--', color='orange', label='Theoretical')
    plt.legend(loc='upper right')
    ax1.grid(True)
    ax1.set_xlabel("Azimuth angle [deg]")
    ax1.set_ylabel("Antenna gain [dBi]")
    plt.title("Meteorological Radar Uniform Antenna Pattern")
    # ax1.set_xlim([0, 40])
    ax1.set_xlim([-90, 90])
    ax1.set_ylim([-60, 0])

    ax2 = fig.add_subplot(122)
    ax2.plot(elevation, gain_elev, 'b--', color='blue', label='Uniform')
    ax2.plot(elevation, theoretical, 'b--', color='orange', label='Theoretical')
    ax2.grid(True)
    ax2.set_xlabel("Elevation angle [deg]")
    ax2.set_ylabel("Antenna gain [dBi]")
    # ax2.set_xlim([0, 40])
    ax2.set_xlim([0, 40])
    ax2.set_ylim([-60, 0])
    plt.title("Meteorological Radar Uniform Antenna Pattern")
    plt.legend(loc='upper right')
    plt.show()
