# -*- coding: utf-8 -*-
"""
Created on Wed May 17 15:47:05 2017

@author: edgar
"""

class ParametersHotspot(object):

    __instance = None

    def __new__(cls):
        """
        This is the Singleton Pattern to ensure that this class will have only
        one instance
        """
        if ParametersHotspot.__instance is None:
            ParametersHotspot.__instance = object.__new__(cls)
        return ParametersHotspot.__instance

    ###########################################################################
    # Number of hotspots per macro cell geographical area
    # Possible values: 1, 2, optional of 4
    num_hotspots_per_cell = 4
    
    ###########################################################################
    # Number of base stations per hotspot
    # Possible values: 4 to 10
    num_stations_per_hotspots = 4
    
    ###########################################################################
    # Ratio of UEs randomly and uniformly dropped within the hotspots
    ue_hotspot_dropping_ratio = 0.67
    
    ###########################################################################
    # Ratio of outdoor UE's in the heterogeneous topology
    ue_outdoor_ratio = 0.8
    
    ###########################################################################
    # Maximum 2D distance between station and hotspot center [m]
    max_dist_station_hotspot = 50
    
    ###########################################################################
    # Maximum 2D distance between UE and hotspot center [m]
    max_dist_ue_hotspot = 70
    
    ###########################################################################
    # Minimum 2D distance between two stations (small cells) in a hotspot [m]
    min_dist_stations = 20
    
    ###########################################################################
    # Minimum 2D distance between station and UE in a hotspot [m]
    min_dist_station_ue = 5
    
    ###########################################################################
    # Minimum 2D distance between macro cell base station and hotspot center [m]
    min_dist_bs_hotspot = 105

    ###########################################################################
    # Minimum 2D distance between macro cell base station and UE [m]
    min_dist_bs_ue = 35

    ###########################################################################
    # Minimum 2D distance between two hotspot centers [m]
    # This is twice the distance between station and hotspot center
    min_dist_hotspots = 2*max_dist_station_hotspot