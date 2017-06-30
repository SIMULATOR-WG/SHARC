# -*- coding: utf-8 -*-
"""
Created on Tue Feb 14 12:48:58 2017

@author: edgar
"""

from abc import ABCMeta, abstractmethod
import numpy as np
 
class Topology(object):
    
    __metaclass__ = ABCMeta
    
    def __init__(self, intersite_distance: float, cell_radius: float):
        self.intersite_distance = intersite_distance
        self.cell_radius = cell_radius
        
        # Coordinates of the base stations. In this context, each base station
        # is equivalent to a sector (hexagon) in the macrocell topology
        self.x = np.empty(0)
        self.y = np.empty(0)
        self.azimuth = np.empty(0)

        # 
        self.calculate_coordinates()
        self.num_base_stations = len(self.x)
    
        
    @abstractmethod
    def calculate_coordinates(self):
        """
        Calculates the coordinates of the stations according to the class
        atributes.
        """        
        pass
        
