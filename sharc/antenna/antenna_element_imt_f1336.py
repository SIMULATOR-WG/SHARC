# -*- coding: utf-8 -*-
"""
Created on Fri Apr 14 14:13:58 2017

@author: Calil
"""

import numpy as np
import sys

from sharc.support.named_tuples import AntennaPar

class AntennaElementImtF1336(object):
    """
    Implements a single element of an IMT antenna array following ITU-R F.1336-4, item 3.1.1
    using parameters from ITU-R M2292

    Attributes
    ----------
        g_max (float): maximum gain of element
        theta_3db (float): vertical 3dB beamwidth of single element [degrees]
        phi_3db (float): horizontal 3dB beamwidth of single element [degrees]
    """

    def __init__(self,par: AntennaPar):
        """
        Constructs an AntennaElementImt object.

        Parameters
        ---------
            param (ParametersAntennaImt): antenna IMT parameters
        """
        self.param = par

        self.g_max = par.element_max_g
        self.phi_3db = par.element_phi_3db
        if par.element_theta_3db > 0:
            self.theta_3db = par.element_theta_3db
        else:
            if self.phi_3db > 120.:
                sys.stderr.write("ERROR\nvertical beamwidth must be givem if horizontal beamwidth > 120 degrees")
                sys.exit(1)
            # calculate based on F1336
            self.theta_3db = (31000 * 10**(-.1 * self.g_max))/self.phi_3db

        # antenna paremeters, according to ITU-R M2292
        self.k_a = .7
        self.k_p = .7
        self.k_h = .7
        self.lambda_k_h = 3 * (1-.5**(-self.k_h))
        self.k_v = .3
        self.incline_factor = 10*np.log10(((180/self.theta_3db)**1.5 * (4**-1.5+self.k_v))/
                                       (1 + 8 * self.k_p)) / np.log10(22.5 / self.theta_3db)
        self.x_k = np.sqrt(1 - .36 * self.k_v)
        self.lambda_k_v = 12 - self.incline_factor * np.log10(4) - 10 * np.log10(4**-1.5 + self.k_v)

        self.g_hr_180 = -12. + 10 * np.log10(1 + 8 * self.k_a) - 15 * np.log10(180/self.theta_3db)
        self.g_hr_0 = 0

    def horizontal_pattern(self, phi: np.array) -> {np.array, float}:
        """
        Calculates the horizontal radiation pattern.

        Parameters
        ----------
            phi (np.array): azimuth angle [degrees]

        Returns
        -------
            a_h (np.array): horizontal radiation pattern gain value
        """
        if type(phi) is not np.ndarray:
            phi_a = np.array([phi])
        else:
            phi_a = phi
        
        x_h = np.abs(phi_a)/self.phi_3db
        gain = np.zeros(np.size(phi_a))

        i0 = np.where(x_h < 0.5)[0]
        gain[i0] = -12 * np.power(x_h[i0], 2)

        i1 = np.where(x_h >= 0.5)[0]
        gain[i1] = -12 * np.power(x_h[i1], 2 - self.k_h) - self.lambda_k_h
        
        gain = np.maximum(gain, self.g_hr_180)
        
        if type(phi) is not np.ndarray:
            gain = gain[0]
        
        return gain

    def vertical_pattern(self, theta: np.array) -> np.array:
        """
        Calculates the vertical radiation pattern.

        Parameters
        ----------
            theta (np.array): elevation angle [degrees]

        Returns
        -------
            a_v (np.array): vertical radiation pattern gain value
        """
        # This correction is needed because the simulator calculates theta 
        # with respect to z-axis and equations of F.1336 assume that theta is
        # calculated with respect to the direction of maximum gain
        if type(theta) is np.ndarray:
            theta_a = theta - 90
        else:
            theta_a = np.array([theta]) - 90
        
        x_v = np.abs(theta_a)/self.theta_3db
        gain = np.zeros(np.size(theta_a))

        i0 = np.where(x_v < self.x_k)[0]
        gain[i0] = -12 * np.power(x_v[i0], 2)
            
        i1 = np.where((x_v >= self.x_k) & (x_v < 4))[0]
        gain[i1] = -12 + 10*np.log10(np.power(x_v[i1], -1.5) + self.k_v)
            
        i2 = np.where((x_v >= 4) & (x_v < 90 / self.theta_3db))[0]
        gain[i2] = - self.lambda_k_v - self.incline_factor * np.log10(x_v[i2])
            
        i3 = np.where(x_v >= (90 / self.theta_3db))[0]
        gain[i3] = self.g_hr_180

        if type(theta) is not np.ndarray:
            gain = gain[0]

        return gain

    def element_pattern(self, phi: np.array, theta: np.array) -> np.array:
        """
        Calculates the element radiation pattern gain.

        Parameters
        ----------
            theta (np.array): elevation angle [degrees]
            phi (np.array): azimuth angle [degrees]

        Returns
        -------
            gain (np.array): element radiation pattern gain value
        """

        gain_hor = self.horizontal_pattern(phi)
        compression_ratio = (gain_hor - self.g_hr_180)/(self.g_hr_0 - self.g_hr_180)
        gain = self.g_max + gain_hor + compression_ratio * self.vertical_pattern(theta)

        return gain

if __name__ == '__main__':

    from sharc.parameters.parameters_antenna_imt import ParametersAntennaImt
    from matplotlib import pyplot as plt

    param = ParametersAntennaImt()

    param.element_max_g = 18
    param.element_phi_3db = 65
    param.element_theta_3db = 0

    antenna = AntennaElementImtF1336( param )

    phi_vec = np.arange(-180, 180, step = 1)
    theta_vec = np.arange(0, 180, step = 1)

    pattern_hor_0deg = antenna.element_pattern(phi_vec, 0)
    pattern_hor_10deg = antenna.element_pattern(phi_vec, 10)
    pattern_hor_30deg = antenna.element_pattern(phi_vec, 30)
    pattern_hor_60deg = antenna.element_pattern(phi_vec, 60)
    pattern_hor_90deg = antenna.element_pattern(phi_vec, 90)
    
    pattern_ver_0deg = antenna.element_pattern(0, theta_vec)
    pattern_ver_30deg = antenna.element_pattern(30, theta_vec)
    pattern_ver_60deg = antenna.element_pattern(60, theta_vec)
    pattern_ver_90deg = antenna.element_pattern(90, theta_vec)
    pattern_ver_120deg = antenna.element_pattern(120, theta_vec)

    plt.figure(1)
    plt.plot(phi_vec, pattern_hor_0deg, label = 'elevation = 0 degrees')
    plt.plot(phi_vec, pattern_hor_10deg, label = 'elevation = 10 degrees')
    plt.plot(phi_vec, pattern_hor_30deg, label = 'elevation = 30 degrees')
    plt.plot(phi_vec, pattern_hor_60deg, label = 'elevation = 60 degrees')
    plt.plot(phi_vec, pattern_hor_90deg, label = 'elevation = 90 degrees')

    plt.title('horizontal pattern')
    plt.xlabel('azimuth (degrees)')
    plt.ylabel('gain (dBi)')
    plt.legend()

    plt.figure(2)
    plt.plot(theta_vec, pattern_ver_0deg, label='azimuth = 0 degrees')
    plt.plot(theta_vec, pattern_ver_30deg, label='azimuth = 30 degrees')
    plt.plot(theta_vec, pattern_ver_60deg, label='azimuth = 60 degrees')
    plt.plot(theta_vec, pattern_ver_90deg, label='azimuth = 90 degrees')
    plt.plot(theta_vec, pattern_ver_120deg, label='azimuth = 120 degrees')

    plt.title('vertical pattern')
    plt.xlabel('elevation (degrees)')
    plt.ylabel('gain (dBi)')
    plt.legend()

    plt.show()
