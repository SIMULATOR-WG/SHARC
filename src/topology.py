# -*- coding: utf-8 -*-
"""
Created on Tue Feb 14 12:48:58 2017

@author: edgar
"""

from abc import ABCMeta, abstractmethod
import numpy as np
import typing
 
class Topology(object):
    
    __metaclass__ = ABCMeta
    
    def __init__(self, intersite_distance: float, cell_radius: float, 
                 num_clusters: int):
        self.__intersite_distance = intersite_distance
        self.__cell_radius = cell_radius
        self.__num_clusters = num_clusters
        self._calculate_coordinates()
        self._calculate_limits()
    
    @abstractmethod
    def _calculate_coordinates(self):
        self.__x_coord = np.empty(0)
        self.__y_coord = np.empty(0)
        
    @abstractmethod
    def _calculate_limits(self):
        self.__x_min = 0
        self.__x_max = 0
        self.__y_min = 0
        self.__y_max = 0
    
    def get_coordinates(self) -> typing.Tuple[np.array, np.array]:
        return (self.__x_coord, self.__y_coord)
        
    def get_x_limits(self) -> typing.Tuple[float, float]:
        return (self.__x_min, self.__x_max)

    def get_y_limits(self) -> typing.Tuple[float, float]:
        return (self.__y_min, self.__y_max)    
        
    def set_coordinates(self, x_coord: np.array, y_coord: np.array):
        self.__x_coord = np.asarray(x_coord)
        self.__y_coord = np.asarray(y_coord)
        
    def set_x_limits(self, x_min: float, x_max: float):
        self.__x_min = x_min
        self.__x_max = x_max

    def set_y_limits(self, y_min: float, y_max: float):
        self.__y_min = y_min
        self.__y_max = y_max         
        
    def set_intersite_distance(self, intersite_distance: float):
        self.__intersite_distance = intersite_distance
        self._calculate_coordinates()
        self._calculate_limits()        
        
    def get_intersite_distance(self) -> float:
        return self.__intersite_distance     
        
    def set_cell_radius(self, cell_radius: float):
        self.__cell_radius = cell_radius
        self._calculate_coordinates()
        self._calculate_limits()        
        
    def get_cell_radius(self) -> float:
        return self.__cell_radius         
        
    def set_num_clusters(self, num_clusters: int):
        self.__num_clusters = num_clusters
        self._calculate_coordinates()
        self._calculate_limits()        
        
    def get_num_clusters(self) -> int:
        return self.__num_clusters        