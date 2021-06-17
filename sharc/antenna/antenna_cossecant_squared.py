# -*- coding: utf-8 -*-
"""
@Created: Luciano Camilo on Thu Marc 25 10:34:00 2021

"""
from sharc.antenna.antenna import Antenna
from sharc.parameters.parameters_arns import ParametersArns

import numpy as np


class AntennaCossecantSquared(Antenna):
    """
    Implements a Cossecant-squared (CSC2) antenna pattern of an Aeronautical Surveillance Radar antenna following
    ITU-R M.1851-1, Table 4 and Figure 12.
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
        self.maximum_csc = param.csc2_angle
        self.highbeam_csc2 = param.highbeam_csc2

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
        #print(f'Angulo theta de incidencia: {theta_l}')
        #print(f'Ganho da antena no angulo: {gain}')
        return gain

    def get_gain_az(self, phi: np.array) -> np.array:
        """
        Returns the antenna gain in the azimuth plane for the given direction.
        """
        const = np.pi * 68.8
        g = np.zeros(len(phi))
        mi = np.zeros(len(phi))

        for i in range(len(phi)):

            if -self.phi_3db <= phi[i] <= self.phi_3db:
                mi[i] = const * np.sin(np.radians(phi[i])) / self.phi_3db
                g[i] = 20 * np.log10((np.pi / 2) * ((np.cos(mi[i])) / ((np.pi / 2) ** 2 - (mi[i]) ** 2))) + 4.32 - 0.392

            if self.phi_3db <= phi[i] <= 180 or -self.phi_3db >= phi[i] >= -180:
                g[i] = -17.51 * np.log(2.33 * (np.abs((phi[i])) / self.phi_3db)) - 4.32  # cos2

                if g[i] < -50:
                    g[i] = -50
        return g

    def get_gain_elev(self, theta: np.array) -> np.array:
        """
        Returns the antenna gain in the elevation plane for the given direction.
        """

        const = np.pi * 50.8
        g = np.zeros(len(theta))
        mi = np.zeros(len(theta))
        theta = theta-self.highbeam_csc2
        self.maximum_csc = self.maximum_csc-self.highbeam_csc2

        for i in range(len(theta)):
            if -self.theta_3db / 0.88 <= theta[i] <= self.theta_3db:
                mi[i] = const * np.sin(np.radians(theta[i])) / np.radians(self.theta_3db)
                g[i] = 20 * np.log10(np.sin(np.radians(mi[i])) / mi[i]) + 35.1623
            if self.theta_3db <= theta[i] <= self.maximum_csc:
                g1 = ((np.sin(np.radians(const * np.sin(np.radians(self.theta_3db)) / np.radians(self.theta_3db))))
                      / (const * np.sin(np.radians(self.theta_3db)) / np.radians(self.theta_3db))) + 0.12
                g[i] = 20 * np.log10(
                    g1 * (((1 / np.sin(np.radians(theta[i]))) / (1 / np.sin(np.radians(self.theta_3db)))) ** 2))
            if self.maximum_csc <= theta[i] <= 90:
                g[i] = -55
            if theta[i] < -self.theta_3db / 0.88 or theta[i] >= 90:
                g[i] = -55

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
    param.beamwidth_el = 4
    param.beamwidth_az = 1.5
    param.csc2_angle = 35
    param.highbeam_csc2 = 0

    azimuth = np.linspace(-90, 90, num=1500)
    elevation = np.linspace(-90, 90, num=1500)

    gain = 0
    azimuth_l = 0
    elevation_l = 0

    antenna = AntennaCossecantSquared(param)

    gain_az = gain + antenna.get_gain_az(azimuth)
    gain_elev = gain + antenna.get_gain_elev(elevation)

    fig = plt.figure(figsize=(15, 5), facecolor='w', edgecolor='k')

    ax1 = fig.add_subplot(121)
    #ax1.plot(azimuth, gain_az, 'b--', color='blue', label='Cosine Antenna Pattern')
    ax1.semilogx(azimuth, gain_az, 'b--', color='blue', label='Cosine Antenna Pattern')
    plt.legend()
    ax1.grid(True)
    ax1.set_xlabel("Azimuth angle [deg]")
    ax1.set_ylabel("Antenna gain [dBi]")
    plt.title("Aeronautical Surveillance Radar Cosine Antenna Pattern")
    ax1.set_xlim([0, 30])
    # ax1.set_ylim([-90, 90])

    ax2 = fig.add_subplot(122)
    ax2.plot(elevation, gain_elev, 'b--', color='darkorange', label='$CSC^2$ Antenna Pattern')
    ax2.grid(True)
    ax2.set_xlabel("Elevation angle [deg]")
    ax2.set_ylabel("Antenna gain [dBi]")
    ax2.set_xlim([-5, 50])
    # ax2.set_ylim([-30, 40])
    plt.title("Aeronautical Surveillance Radar CSC2 Antenna Pattern")
    plt.legend()
    plt.show()
