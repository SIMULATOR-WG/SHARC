#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 10 15:19:48 2017

@author: carlosrodriguez
"""
import math

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
        diameter (float): diamteer of earth station antenna [m]
        frequency (float): frequency of operation [MHz]
    """
    
    def __init__(self, diameter: float, frequency: float):
        """
        Constructs an AntennaImt object.
        
        Parameters
        ---------
        diameter: diameter of the earth station antenna
        frequency: operation frequency of the antena of the earth station
        """
        self.__diameter = diameter
        self.__frequency = frequency
        
    def get_gain(self, phi: float, theta: float) -> float:
        wavelength = 3e8 / (self.__frequency * 1000000)
        d_to_wavel = self.__diameter/wavelength
        phimin1 = 15.85 * math.pow(d_to_wavel, -0.6)
        phimin2 = 118 * math.pow(d_to_wavel, -1.06)
        if phimin1 > phimin2:
            phimin = phimin1
        else:
            phimin = phimin2
                
        
        if d_to_wavel >= 46.8:
            if phi >= phimin and phi <= 7:
                return 29 + 3 * math.pow(math.sin(theta * math.pi / 180) , 2) - 25 * math.log10(phi)
            elif phi > 7 and phi <= 9.2:
                return 7.9 + (3 * math.pow(math.sin(theta * math.pi / 180),2)) * (9.2 - phi) / 2.2
            elif phi > 9.2 and phi <= 48:
                return 32 - 25 * math.log10(phi)
            else:
                return -10  
        elif d_to_wavel < 46.8 and d_to_wavel >= 15:
            if phi >= phimin and phi <= 7:
                return 29 + 3 * math.pow(math.sin(theta * math.pi / 180),2) - 25 * math.log10(phi)
            elif phi > 7 and phi <= 9.2: 
                return 7.9 + (3 * math.pow(math.sin(theta * math.pi / 180)),2) * (9.2 - phi) / 2.2
            elif phi > 9.2 and phi <= 30.2:
                return 32 - 25 * math.log10(phi)
            elif phi > 30.2 and phi <= 70:
                return -5
            else:
                return 0             
        else:
            return 0
                