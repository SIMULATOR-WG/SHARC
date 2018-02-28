# -*- coding: utf-8 -*-
"""
Created on Fri Apr 14 14:13:58 2017

@author: Calil
"""

import numpy as np
import sys

from sharc.support.named_tuples import AntennaPar
from sharc.antenna.antenna import Antenna

class AntennaSectorF1336(Antenna):
    """
    Implements a single element of an IMT antenna array following ITU-R F.1336-4, item 3.1.1
    using parameters from ITU-R M2292

    Attributes
    ----------
        g_max (float): maximum gain of element
        theta_3db (float): vertical 3dB beamwidth of single element [degrees]
        phi_3db (float): horizontal 3dB beamwidth of single element [degrees]
        am (float): front-to-back ratio
        sla_v (float): element vertical sidelobe attenuation
    """

    def __init__(self,par: AntennaPar,downtilt_deg: float, elevation:float, azimuth: float):
        """
        Constructs an AntennaElementImt object.

        Parameters
        ---------
            param (ParametersAntennaImt): antenna IMT parameters
        """
        super().__init__()
        self.param = par
        
        self.elevation = elevation
        self.azimuth = azimuth

        self.g_max = par.element_max_g
        self.downtilt_rad = downtilt_deg / 180 * np.pi
        self.phi_deg_3db = par.element_phi_3db
        if par.element_theta_3db > 0:
            self.theta_deg_3db = par.element_theta_3db
        else:
            if self.phi_deg_3db > 120.:
                sys.stderr.write("ERROR\nvertical beamwidth must be givem if horizontal beamwidth > 120 degrees")
                sys.exit(1)
            # calculate based on F1336
            self.theta_deg_3db = (31000 * 10**(-.1 * self.g_max))/self.phi_deg_3db

        # antenna paremeters, according to ITU-R M2292
        self.k_a = .7
        self.k_p = .7
        self.k_h = .7
        self.lambda_k_h = 3 * (1-.5**(-self.k_h))
        self.k_v = .3
        self.incline_factor = np.log10(((180/self.theta_deg_3db)**1.5 * (4**-1.5+self.k_v))/
                                       (1 + 8 * self.k_p)) / np.log10(22.5 / self.theta_deg_3db)
        self.x_k = np.sqrt(1 - .36 * self.k_v)
        self.lambda_k_v = 12 - self.incline_factor * np.log10(4) - 10 * np.log10(4**-1.5 + self.k_v)

        self.g_hr_180 = -12. + 10 * np.log10(1 + 8 * self.k_a) - 15 * np.log10(180/self.theta_deg_3db)
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
        x_h = abs(phi)/self.phi_deg_3db
        if x_h < 0.5:
            gain = -12 * x_h ** 2
        else:
            gain = -12 * x_h ** (2 - self.k_h) - self.lambda_k_h
        gain = max(gain, self.g_hr_180)

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
        x_v = abs(theta)/self.theta_deg_3db

        if x_v < self.x_k:
            gain = -12 * x_v ** 2
        elif x_v < 4:
            gain = -12 + 10*np.log10(x_v**-1.5 + self.k_v)
        elif x_v < 90 / self.theta_deg_3db:
            gain = - self.lambda_k_v - self.incline_factor * np.log10(x_v)
        else:
            gain = self.g_hr_180

        return gain

    def calculate_gain(self, *args, **kwargs) -> np.array:
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
        theta_v = np.asarray(kwargs["theta_vec"])
        phi_v = np.asarray(kwargs["phi_vec"])
        
        phi, theta = self.to_local_coord(phi_v,theta_v)

        # recalculate angles considering mechanical tilt (eqs 3b/3c)
        theta_rad = -theta / 180 * np.pi
        phi_rad = phi / 180 * np.pi
        new_theta_rad = np.arcsin(np.sin(theta_rad) * np.cos(self.downtilt_rad) +
                                  np.cos(theta_rad) * np.cos(phi_rad) * np.sin(self.downtilt_rad))
        cos = (-np.sin(theta_rad) * np.sin(self.downtilt_rad) +
                np.cos(theta_rad) * np.cos(phi_rad) * np.cos(self.downtilt_rad))/np.cos(new_theta_rad)

        # to avoid numerical errors, as sometimes cosines are slightly out of bounds
        if cos > 1:
            cos = 1
        elif cos < -1:
            cos = -1

        phi_rad = np.arccos(cos)
        theta = new_theta_rad / np.pi * 180
        phi = phi_rad / np.pi * 180

        #theta = theta - self.downtilt_rad * 180 / np.pi
        gain_hor = self.horizontal_pattern(phi)
        compression_ratio = (gain_hor - self.g_hr_180)/(self.g_hr_0 - self.g_hr_180)
        gain = self.g_max + gain_hor + compression_ratio * self.vertical_pattern(theta)

        return gain
    
    def to_local_coord(self,phi,theta):
        """
        Receives theta with reference to z axis, converts it to reference in
        x axis and converts to local coordinate system
        """
        lo_theta = np.ravel(np.array([theta + self.elevation]))
        lo_phi = np.ravel(np.array([phi - self.azimuth]))
        
        lo_theta = 90 - np.ravel(np.mod(np.array([lo_theta]),360))
        
        ofb_theta = np.where(np.logical_or(lo_theta < -90,lo_theta > 90))
        lo_theta[ofb_theta] = np.sign(lo_theta[ofb_theta])*180 - lo_theta[ofb_theta]
        lo_phi[ofb_theta] = lo_phi[ofb_theta] + 180
        
        ofb_phi = np.where(np.logical_or(lo_phi < -180,lo_phi > 180))
        lo_phi[ofb_phi] = np.mod(lo_phi[ofb_phi],360)
        ofb_phi = np.where(lo_phi > 180)
        lo_phi[ofb_phi] = lo_phi[ofb_phi] - 360
        
        return lo_phi, lo_theta

if __name__ == '__main__':

    from sharc.parameters.parameters_antenna_imt import ParametersAntennaImt
    from matplotlib import pyplot as plt

    param = ParametersAntennaImt()

    param.element_max_g = 15
    param.element_phi_3db = 65
    param.element_theta_3db = 0

    # 0 degrees tilt
    elevation = 0
    azimuth = 0
    downtilt_deg = 0

    antenna = AntennaSectorF1336(param,downtilt_deg,elevation,azimuth)

    phi_v = np.arange(-180,180, step = 5)
    theta_v = np.arange(0,180, step = 3)

    pattern_hor_0deg = np.zeros(phi_v.shape)
    pattern_hor_10deg = np.zeros(phi_v.shape)
    pattern_hor_30deg = np.zeros(phi_v.shape)
    pattern_hor_60deg = np.zeros(phi_v.shape)

    pattern_ver_0deg = np.zeros(theta_v.shape)
    pattern_ver_30deg = np.zeros(theta_v.shape)
    pattern_ver_60deg = np.zeros(theta_v.shape)
    pattern_ver_90deg = np.zeros(theta_v.shape)
    pattern_ver_120deg = np.zeros(theta_v.shape)

    for phi, index in zip(phi_v, range(len(phi_v))):
        pattern_hor_0deg[index] = antenna.calculate_gain(phi_vec=phi,theta_vec=0)
        pattern_hor_10deg[index] = antenna.calculate_gain(phi_vec=phi,theta_vec=10)
        pattern_hor_30deg[index] = antenna.calculate_gain(phi_vec=phi,theta_vec=30)
        pattern_hor_60deg[index] = antenna.calculate_gain(phi_vec=phi,theta_vec=60)

    plt.figure(1)
    plt.plot(phi_v, pattern_hor_0deg, label = 'elevation = 0 degrees')
    plt.plot(phi_v, pattern_hor_10deg, label = 'elevation = 10 degrees')
    plt.plot(phi_v, pattern_hor_30deg, label = 'elevation = 30 degrees')
    plt.plot(phi_v, pattern_hor_60deg, label = 'elevation = 60 degrees')

    plt.title('downtilt = 0 degrees')
    plt.xlabel ('azimuth (degrees)')
    plt.ylabel ('gain (dBi)')

    plt.legend()

    for theta, index in zip(theta_v, range(len(theta_v))):
        pattern_ver_0deg[index] = antenna.calculate_gain(phi_vec=0,theta_vec=theta)
        pattern_ver_30deg[index] = antenna.calculate_gain(phi_vec=30,theta_vec=theta)
        pattern_ver_60deg[index] = antenna.calculate_gain(phi_vec=60,theta_vec=theta)
        pattern_ver_90deg[index] = antenna.calculate_gain(phi_vec=90,theta_vec=theta)
        pattern_ver_120deg[index] = antenna.calculate_gain(phi_vec=120,theta_vec=theta)

    plt.figure(2)
    plt.plot(theta_v, pattern_ver_0deg, label='azimuth = 0 degrees')
    plt.plot(theta_v, pattern_ver_30deg, label='azimuth = 30 degrees')
    plt.plot(theta_v, pattern_ver_60deg, label='azimuth = 60 degrees')
    plt.plot(theta_v, pattern_ver_90deg, label='azimuth = 90 degrees')
    plt.plot(theta_v, pattern_ver_120deg, label='azimuth = 120 degrees')

    plt.title('downtilt = 0 degrees')
    plt.xlabel('elevation (degrees)')
    plt.ylabel('gain (dBi)')

    plt.legend()

    # x degrees tilt
    downtilt_deg = -10
    antenna = AntennaSectorF1336(param,downtilt_deg,elevation,azimuth)

    for phi, index in zip(phi_v, range(len(phi_v))):
        pattern_hor_0deg[index] = antenna.calculate_gain(phi_vec=phi,theta_vec=0)
        pattern_hor_10deg[index] = antenna.calculate_gain(phi_vec=phi,theta_vec=10)
        pattern_hor_30deg[index] = antenna.calculate_gain(phi_vec=phi,theta_vec=30)
        pattern_hor_60deg[index] = antenna.calculate_gain(phi_vec=phi,theta_vec=60)

    plt.figure(3)
    plt.plot(phi_v, pattern_hor_0deg, label='0 degrees')
    plt.plot(phi_v, pattern_hor_10deg, label='10 degrees')
    plt.plot(phi_v, pattern_hor_30deg, label='30 degrees')
    plt.plot(phi_v, pattern_hor_60deg, label='60 degrees')

    plt.title('downtilt = {} degrees'.format(downtilt_deg))
    plt.xlabel ('azimuth (degrees)')
    plt.ylabel ('gain (dBi)')
    plt.legend()

    for theta, index in zip(theta_v, range(len(theta_v))):
        pattern_ver_0deg[index] = antenna.calculate_gain(phi_vec=0,theta_vec=theta)
        pattern_ver_30deg[index] = antenna.calculate_gain(phi_vec=30,theta_vec=theta)
        pattern_ver_60deg[index] = antenna.calculate_gain(phi_vec=60,theta_vec=theta)
        pattern_ver_90deg[index] = antenna.calculate_gain(phi_vec=90,theta_vec=theta)
        pattern_ver_120deg[index] = antenna.calculate_gain(phi_vec=120,theta_vec=theta)

    plt.figure(4)
    plt.plot(theta_v, pattern_ver_0deg, label='azimuth = 0 degrees')
    plt.plot(theta_v, pattern_ver_30deg, label='azimuth = 30 degrees')
    plt.plot(theta_v, pattern_ver_60deg, label='azimuth = 60 degrees')
    plt.plot(theta_v, pattern_ver_90deg, label='azimuth = 90 degrees')
    plt.plot(theta_v, pattern_ver_120deg, label='azimuth = 120 degrees')

    plt.title('downtilt = {} degrees'.format(downtilt_deg))
    plt.xlabel('elevation (degrees)')
    plt.ylabel('gain (dBi)')

    plt.legend()
    plt.show()
