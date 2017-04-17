# -*- coding: utf-8 -*-
"""
Created on Thu Apr 13 13:16:02 2017

@author: edgar
"""

class ParametersFss(object):

    __instance = None
    
    def __new__(cls):
        """
        This is the Singleton Pattern to ensure that this class will have only
        one instance
        """
        if ParametersFss.__instance is None:
            ParametersFss.__instance = object.__new__(cls)
        return ParametersFss.__instance    
    
    ###########################################################################
    # satellite center frequency [MHz]
    sat_frequency = 27250

    ###########################################################################
    # satellite bandwidth [MHz]
    sat_bandwidth = 200
        
    ###########################################################################
    # satellite altitude [m]
    sat_height = 35786000
    
    ###########################################################################
    # Satellite peak reeive antenna gain [dBi]
    sat_rx_antenna_gain = 51
    
    ###########################################################################
    # System receive noise temperature [K]
    sat_noise_temperature = 950

    ###########################################################################
    # Interference protection criteria [dB]
    sat_interference_noise_ratio = -12.2

    ###########################################################################
    # Satellite antenna pattern in the fixed-satellite service
    sat_antenna_pattern = "ITU-R S.672-4"
    
    ###########################################################################
    # 
    BOLTZMANN_CONSTANT = 1.38064852e-23