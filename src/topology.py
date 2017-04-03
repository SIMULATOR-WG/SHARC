# -*- coding: utf-8 -*-
"""
Created on Tue Feb 14 12:48:58 2017

@author: edgar
"""

from abc import ABCMeta, abstractmethod
import numpy as np
 
class Topology(object):
    
    __metaclass__ = ABCMeta
    
    def __init__(self, intersite_distance: float, cell_radius: float, 
                 num_clusters: int, allowed_num_clusters: list):
        self.allowed_num_clusters = allowed_num_clusters
        self._intersite_distance = intersite_distance
        self._cell_radius = cell_radius
        self._num_clusters = num_clusters
        self.x = np.empty(0)
        self.y = np.empty(0)        
        self.x_min = 0
        self.x_max = 0
        self.y_min = 0
        self.y_max = 0
        self._calculate_coordinates()
        self._calculate_limits()
    
    @property
    def intersite_distance(self):
        return self._intersite_distance
        
    @intersite_distance.setter
    def intersite_distance(self, value):
        """
        Sets class atribute and recalculates coordinates and limits
        """
        self._intersite_distance = value
        self._calculate_coordinates()
        self._calculate_limits()        

    @property
    def cell_radius(self):
        return self._cell_radius

    @cell_radius.setter
    def cell_radius(self, value):
        """
        Sets class atribute and recalculates coordinates and limits
        """
        self._cell_radius = value
        self._calculate_coordinates()
        self._calculate_limits() 
        
    @property
    def num_clusters(self):
        return self._num_clusters        
    
    @num_clusters.setter
    def num_clusters(self, value):
        """
        Sets class atribute and recalculates coordinates and limits
        """
        if value not in self.allowed_num_clusters:
            error_message = "invalid number of clusters ({})".format(value)
            raise ValueError(error_message) 
        self._num_clusters = value
        self._calculate_coordinates()
        self._calculate_limits()   
        
    @abstractmethod
    def _calculate_coordinates(self):
        """
        Calculates the coordinates of the stations according to the class
        atributes.
        """        
        pass
        
    def _calculate_limits(self):
        """
        Calculates the coordinates of the scenario's borders
        """        
        self.x_min = np.min(self.x) - self.cell_radius
        self.x_max = np.max(self.x) + self.cell_radius
        self.y_min = np.min(self.y) - self.cell_radius
        self.y_max = np.max(self.y) + self.cell_radius        
        pass