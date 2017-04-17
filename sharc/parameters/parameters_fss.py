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
    # satellite altitude [m] and latitude [deg]
    sat_altitude = 35786000
    sat_lat_deg = 0

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
    # IMT parameters relevant to the satellite system
    imt_altitude = 1000   # altitude of IMT system (in meters)
    imt_lat_deg = -15.795038   # latitude of IMT system (in degrees)
    imt_long_diff_deg = 3 # difference between longitudes of IMT and satellite system
                          # positive if space-station is to the East of earth-station


    ###########################################################################
    # Channel parameters
    # channel model, possible values are "FSPL" (free-space path loss),
    #                                    "P619" (ITU-R P.619-1)
    channel_model = "P619"
    ###########################################################################

    ###########################################################################
    # Constants
    BOLTZMANN_CONSTANT = 1.38064852e-23
    EARTH_RADIUS = 6371000
