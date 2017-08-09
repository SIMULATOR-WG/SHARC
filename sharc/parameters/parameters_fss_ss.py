# -*- coding: utf-8 -*-
"""
Created on Thu Apr 13 13:16:02 2017

@author: edgar
"""

class ParametersFssSs(object):

    __instance = None

    def __new__(cls):
        """
        This is the Singleton Pattern to ensure that this class will have only
        one instance
        """
        if ParametersFssSs.__instance is None:
            ParametersFssSs.__instance = object.__new__(cls)
        return ParametersFssSs.__instance

    ###########################################################################
    # satellite center frequency [MHz]
    frequency = 27250

    ###########################################################################
    # satellite bandwidth [MHz]
    bandwidth = 200

    ###########################################################################
    # satellite altitude [m] and latitude [deg]
    altitude = 35780000
    lat_deg = 0

    ###########################################################################
    # Elevation angle [deg]
    elevation = 270
    
    ###########################################################################
    # Azimuth angle [deg]
    azimuth = 0    
    
    ###########################################################################
    # System receive noise temperature [K]
    noise_temperature = 950

    ###########################################################################
    # Interference protection criteria [dB]
    interference_noise_ratio = -12.2

    ###########################################################################
    # INR scaling factor (to estimate INR for larger number of interfering stations)
    inr_scaling = 14.95
    
    ###########################################################################
    # Satellite peak reeive antenna gain [dBi]
    antenna_gain = 51
    
    ###########################################################################
    # Satellite antenna pattern in the fixed-satellite service
    # Possible values: "ITU-R S.672-4", "FSS_SS", "OMNI"
    antenna_pattern = "FSS_SS"

    # IMT parameters relevant to the satellite system
    imt_altitude = 0   # altitude of IMT system (in meters)
    imt_lat_deg = 0   # latitude of IMT system (in degrees)
    imt_long_diff_deg = 0 # difference between longitudes of IMT and satellite system
                             # positive if space-station is to the East of earth-station

#    # IMT parameters relevant to the satellite system
#    imt_altitude = 1000   # altitude of IMT system (in meters)
#    imt_lat_deg = -23.5629739   # latitude of IMT system (in degrees)
#    imt_long_diff_deg = (-46.6555132-75) # difference between longitudes of IMT and satellite system
#                             # positive if space-station is to the East of earth-station

    ###########################################################################
    # Channel parameters
    # channel model, possible values are "FSPL" (free-space path loss),
    #                                    "SatelliteSimple" (FSPL + 4 or 24dB (LOS or NLOS)
    #                                    "P619" (ITU-R P.619-1)
    channel_model = "SatelliteSimple"
    line_of_sight_prob = 1 # probability of line-of-sight between UE and satellite

    surf_water_vapour_density = 7.5 #g/m^3
    specific_gaseous_att = 0.1 #db/km
    time_ratio = 0.5 #transmission loss not exceeded for time_ratio*100 % of time

    ###########################################################################

    # The required near-in-side-lobe level (dB) relative to peak gain
    # according to ITU-R S.672-4
    antenna_l_s = -20    
    
    ###########################################################################
    # 3 dB beamwidth angle (3 dB below maximum gain) [degrees]
    antenna_3_dB = 0.65
    
    ###########################################################################
    # Constants
    BOLTZMANN_CONSTANT = 1.38064852e-23
    EARTH_RADIUS = 6371000