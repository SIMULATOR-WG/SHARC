# -*- coding: utf-8 -*-
"""
Created on Mon Nov 27 08:52:28 2017

@author: Calil
"""

from area import area
import numpy as np

class FootprintArea(object):
    """
    Defines a satellite footprint region and calculates its area
    
    
    """
    def __init__(self,bore_lat_deg: float, bore_subsat_long_deg: float, beam:float):
        # Initialize attributes
        self.bore_lat_deg = bore_lat_deg
        self.bore_subsat_long_deg = bore_subsat_long_deg
        self.beam_width = beam
        
        # Calculate tilt
        self.beta = np.arccos(np.cos(np.deg2rad(self.bore_lat_deg))*\
                              np.cos(np.deg2rad(self.bore_subsat_long_deg)))
        self.bore_tilt = np.arctan(np.sin(self.beta)/(6.6235 - np.cos(self.beta)))
        
        
        