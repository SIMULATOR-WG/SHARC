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
    # Number of hotspots per macro cell (sector)
    num_hotspots_per_cell = 1
    
    ###########################################################################
    # Maximum 2D distance between hotspot and UE [m]
    # This is the hotspot radius
    max_dist_hotspot_ue = 100
    
    ###########################################################################
    # Minimum 2D distance between macro cell base station and hotspot [m]
    min_dist_bs_hotspot = 0

    ###########################################################################
    # Minimum 2D distance between two hotspots in the same cell [m]
    # This is twice the distance between hotspot and UE
    min_dist_hotspots = 2*max_dist_hotspot_ue