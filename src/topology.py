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
                 num_clusters: int):
        self._intersite_distance = intersite_distance
        self._cell_radius = cell_radius
        self._num_clusters = num_clusters
        self.__x = np.empty(0)
        self.__y = np.empty(0)        
        self.__x_min = 0
        self.__x_max = 0
        self.__y_min = 0
        self.__y_max = 0
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
        self._num_clusters = value
        self._calculate_coordinates()
        self._calculate_limits()   
        
    @property
    def x(self):
        return self.__x
        
    @x.setter
    def x(self, value):
        self.__x = value
        
    @property
    def y(self):
        return self.__y  
        
    @y.setter
    def y(self, value):
        self.__y = value        
        
    @property
    def x_min(self):
        return self.__x_min

    @x_min.setter
    def x_min(self, value):
        self.__x_min = value
        
    @property
    def x_max(self):
        return self.__x_max

    @x_max.setter
    def x_max(self, value):
        self.__x_max = value
        
    @property
    def y_min(self):
        return self.__y_min

    @y_min.setter
    def y_min(self, value):
        self.__y_min = value
        
    @property
    def y_max(self):
        return self.__y_max        
        
    @y_max.setter
    def y_max(self, value):
        self.__y_max = value
        
    @abstractmethod
    def _calculate_coordinates(self):
        """
        Calculates the coordinates of the stations according to the class
        atributes.
        """        
        pass
        
    @abstractmethod
    def _calculate_limits(self):
        """
        Calculates the coordinates of the scenario's borders
        """        
        pass