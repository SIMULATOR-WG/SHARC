#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 10 15:19:48 2017

@author: carlosrodriguez
"""
import math
import numpy as np

from sharc.antenna.antenna import Antenna

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
    
    def __init__(self, params):
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
        self.azimuth = params.azimuth
        self.elevation = params.elevation
        
    def calculate_gain(self, *args, **kwargs) -> np.array:        
        """
        Calculates the gain of the antenna af arrays of angles.
        
        Parameters
        ----------
            phi_vec (numpy.array): list of azimuth angle [degrees]
            theta_vec (numpy.array) : list of elevation angle [degrees]
            
        Returns
        -------
            gain (numpy.array): gain array in given directions
        """
        
        phi_list = kwargs["phi_vec"]
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
            phi (float): azimuth angle [degrees]
            theta (float) : elevation angle [degrees]
            
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