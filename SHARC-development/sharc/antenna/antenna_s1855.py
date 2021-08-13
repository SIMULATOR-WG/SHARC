#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 10 15:19:48 2017

@author: carlosrodriguez
"""
import math
import numpy as np

from sharc.antenna.antenna import Antenna
from sharc.parameters.parameters_fss_ss import ParametersFssSs

class AntennaS1855(Antenna):
    """"
    Implements amntenna radiation pattern for Earth Station according
    Reccomendation  ITU-R S.1855 (01/2010)
    "Alternative reference radiation pattern for earth station antennas used with
    satellites in the geostationary-satellite orbit for use in coordination and/or
    interference assesment in the frequency range from 2 to 31 GHz"

    Attributes
    ----------
        gain (float): calculated antenna gain in given direction
        diameter (float): diameter of earth station antenna [m]
        frequency (float): frequency of operation [MHz]
    """

    def __init__(self, params: ParametersFssSs):
        """
        Constructs an AntennaS1855 object.

        Parameters
        ---------
            diameter: diameter of the earth station antenna
            frequency: operation frequency of the antena of the earth station
            antenna_gain: peak gain of the antena on the direction of the GSO
        """
        self.diameter = params.diameter
        self.frequency = params.frequency
        self.antenna_gain = params.antenna_gain


    def calculate_gain(self, *args, **kwargs) -> np.array:
        """
        Calculates the gain of the antenna af arrays of angles.

        Parameters
        ----------
            phi_vec (numpy.array): the off-axis angle between the direction of
                interest and the boresight axis [degrees]
            theta_vec (numpy.array) : the angle between the plane containing
                the boresight and the dimension D_GSO, and the plane of interest,
                where the plane of interest passes through the boresight and the
                direction of interest [degrees]

        Returns
        -------
            gain (numpy.array): gain array in given directions
        """

        phi_list = kwargs["off_axis_angle_vec"]
        theta_list = kwargs["theta_vec"]

        gain = np.empty(phi_list.shape, dtype = np.float)

        for i in range(len(phi_list)):
            gain[i] = self.get_gain_pair(phi_list[i], theta_list[i])

        return gain


    def get_gain_pair(self, phi: np.float, theta: np.float) -> np.float:
        """
        Calculates the gain of the antenna of a pair of angles.

        Parameters
        ----------
            phi_vec (numpy.array): the off-axis angle between the direction of
                interest and the boresight axis [degrees]
            theta_vec (numpy.array) : the angle between the plane containing
                the boresight and the dimension D_GSO, and the plane of interest,
                where the plane of interest passes through the boresight and the
                direction of interest [degrees]
        Returns
        -------
            gain (float): gain value in given direction
        """
        gain = None
        wavelength = 3e8 / (self.frequency * 1000000)
        d_to_wavel = self.diameter/wavelength
        phimin1 = 15.85 * math.pow(d_to_wavel, -0.6)
        phimin2 = 118 * math.pow(d_to_wavel, -1.06)
        if phimin1 > phimin2:
            phimin = phimin1
        else:
            phimin = phimin2


        if d_to_wavel >= 46.8:
            if   phi < phimin:
                gain = self.antenna_gain
            elif phi >= phimin and phi <= 7:
                gain = 29 + 3 * np.power(np.sin(theta * math.pi / 180) , 2) - 25 * np.log10(phi)
            elif phi > 7 and phi <= 9.2:
                gain = 7.9 + (3 * np.power(np.sin(theta * np.pi / 180),2)) * (9.2 - phi) / 2.2
            elif phi > 9.2 and phi <= 48:
                gain = 32 - 25 * np.log10(phi)
            else:
                return -10
        elif d_to_wavel < 46.8 and d_to_wavel >= 15:
            if   phi < phimin:
                gain = self.antenna_gain
            elif phi >= phimin and phi <= 7:
                gain = 29 + 3 * np.pow(np.sin(theta * np.pi / 180),2) - 25 * np.log10(phi)
            elif phi > 7 and phi <= 9.2:
                gain = 7.9 + (3 * np.pow(np.sin(theta * np.pi / 180)),2) * (9.2 - phi) / 2.2
            elif phi > 9.2 and phi <= 30.2:
                gain = 32 - 25 * np.log10(phi)
            elif phi > 30.2 and phi <= 70:
                gain = -5
            else:
                gain = 0
        else:
            gain = 0
        return gain


if __name__ == '__main__':

    import matplotlib.pyplot as plt

    params_fss_ss = ParametersFssSs()
    params_fss_ss.diameter = 9.1
    params_fss_ss.frequency = 24250
    params_fss_ss.antenna_gain = 62

    antenna = AntennaS1855(params_fss_ss)

    # Plot radiation pattern for theta = 90 degrees
    off_axis_angle_vec = np.linspace(0.01, 180, num = 10000)
    theta_90 = 90*np.ones(off_axis_angle_vec.shape)
    theta = 0*np.ones(off_axis_angle_vec.shape)

    gain_90 = antenna.calculate_gain(off_axis_angle_vec = off_axis_angle_vec, theta_vec = theta_90)
    gain = antenna.calculate_gain(off_axis_angle_vec = off_axis_angle_vec, theta_vec = theta)

    fig = plt.figure(figsize=(8,7), facecolor='w', edgecolor='k')
    plt.semilogx(off_axis_angle_vec, gain_90, "-b", label = "$\\theta = 90$ deg")
    plt.semilogx(off_axis_angle_vec, gain, "-r", label = "$\\theta = 0$ deg")

    plt.xlabel("Off-axis angle, $\\varphi$ [deg]")
    plt.ylabel("Gain [dBi]")
    plt.title("ITU-R S.1855 Radiation Pattern")
    plt.ylim((-15, 65))
    plt.xlim((0.01, 180))
    plt.legend(loc="upper right")

    plt.grid()
    plt.show()

#    gain_dir = "/Users/edgar/Desktop/"
#    file_name = gain_dir + "S1855.png"
#    plt.savefig(file_name)
