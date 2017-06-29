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
        self.x_min = 0
        self.x_max = 0
        self.y_min = 0
        self.y_max = 0
        
        # 
        self.calculate_coordinates()
        self.calculate_limits()
    
        
    @abstractmethod
    def calculate_coordinates(self):
        """
        Calculates the coordinates of the stations according to the class
        atributes.
        """        
        pass
        
    def calculate_limits(self):
        """
        Calculates the coordinates of the scenario's borders
        """        
        self.x_min = np.min(self.x) - self.cell_radius
        self.x_max = np.max(self.x) + self.cell_radius
        self.y_min = np.min(self.y) - self.cell_radius
        self.y_max = np.max(self.y) + self.cell_radius        
        pass