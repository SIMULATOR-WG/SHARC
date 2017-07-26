# -*- coding: utf-8 -*-
"""
Created on Mon Jul 24 15:30:49 2017

@author: edgar
"""

class ParametersFssEs(object):
    """
    Simulation parameters for FSS Earth Station.
    """
    
    __instance = None

    def __new__(cls):
        """
        This is the Singleton Pattern to ensure that this class will have only
        one instance
        """
        if ParametersFssEs.__instance is None:
            ParametersFssEs.__instance = object.__new__(cls)
        return ParametersFssEs.__instance

    ###########################################################################
    # x-y coordinates [m]
    x = -5000
    y = 0
        
    ###########################################################################
    # antenna height [m]
    height = 10
    
    ###########################################################################
    # Elevation angle [deg]
    elevation = 20
    
    ###########################################################################
    # Azimuth angle [deg]
    azimuth = 0

    ###########################################################################
    # center frequency [MHz]
    frequency = 27250

    ###########################################################################
    # bandwidth [MHz]
    bandwidth = 200

    ###########################################################################
    # Peak transmit power spectral density (clear sky) [dBW/Hz]
    tx_power_density = -68.3
    
    ###########################################################################
    # antenna peak gain [dBi]
    antenna_gain = 62.8
    
    ###########################################################################
    # Antenna pattern 
    antenna_pattern = "ITU-R S.1855"

    ###########################################################################
    # Channel parameters
    # channel model, possible values are "FSPL" (free-space path loss),
    #                                    "SatelliteSimple" (FSPL + 4 or 24dB (LOS or NLOS)
    #                                    "P619" (ITU-R P.619-1)
    channel_model = "FSPL"

    ###########################################################################
    # Line of sight probability between FSS and IMT stations
    line_of_sight_prob = 1 

    ###########################################################################
    # Constants
    BOLTZMANN_CONSTANT = 1.38064852e-23
    EARTH_RADIUS = 6371000