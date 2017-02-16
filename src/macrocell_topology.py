# -*- coding: utf-8 -*-
"""
Created on Tue Feb 14 12:51:22 2017

@author: edgar
"""

from topology import Topology

import math
import numpy as np
import typing

class MacrocellTopology(Topology):
    """
    Generates the coordinates of the stations based on the macrocell network
    topology.
    """
    
    def __init__(self, intersite_distance: float, num_clusters: int):
        """
        Constructor method that sets the parameters and already calls the 
        calculation methods.
        
        Parameters
        ----------
            intersite_distance : Distance between stations
            num_clusters : Number of cluters, should be 1 or 7
        """
        self.__ALLOWED_NUM_CLUSTERS = [1,7]
        self.set_intersite_distance(intersite_distance)
        self.set_num_clusters(num_clusters)
        self.__calculate_coordinates()
        self.__calculate_limits()

    def __calculate_coordinates(self):
        """
        Calculates the coordinates of the stations according to the inter-site
        distance parameter.
        """
        d = self.get_intersite_distance()
        h = (d/3)*math.sqrt(3)/2
        # these are the coordinates of the central cluster
        x_central = np.array([0, d, d/2, -d/2, -d, -d/2, 
                         d/2, 2*d, 3*d/2, d, 0, -d, 
                         -3*d/2, -2*d, -3*d/2, -d, 0, d, 3*d/2])
        y_central = np.array([0, 0, 3*h, 3*h, 0, -3*h, 
                         -3*h, 0, 3*h, 6*h, 6*h, 6*h, 
                         3*h, 0, -3*h, -6*h, -6*h, -6*h, -3*h])
        self.__x_coord = np.copy(x_central)
        self.__y_coord = np.copy(y_central)
        # other clusters are calculated by shifting the central cluster
        if self.get_num_clusters() == 7:
            x_shift = np.array([7*d/2, -d/2, -4*d, -7*d/2, d/2, 4*d])
            y_shift = np.array([9*h, 15*h, 6*h, -9*h, -15*h, -6*h])
            for x, y in zip(x_shift, y_shift):
                self.__x_coord = np.concatenate((self.__x_coord, x_central+x))
                self.__y_coord = np.concatenate((self.__y_coord, y_central+y))        
    
    def __calculate_limits(self):
        """
        Calculates the coordinates of the scenario's borders
        """
        cell_radius = self.get_intersite_distance()*2/3
        self.__x_min = np.min(self.__x_coord) - cell_radius
        self.__x_max = np.max(self.__x_coord) + cell_radius
        self.__y_min = np.min(self.__y_coord) - cell_radius
        self.__y_max = np.max(self.__y_coord) + cell_radius
        
    def set_intersite_distance(self, intersite_distance: float):
        self.__intersite_distance = intersite_distance
        
    def get_intersite_distance(self) -> float:
        return self.__intersite_distance       
        
    def set_num_clusters(self, num_clusters: int):
        """
        TODO: catch the exception in model and display it in console/log
        """
        if num_clusters not in self.__ALLOWED_NUM_CLUSTERS:
            error_message = "invalid number of clusters ({})".format(num_clusters)
            raise ValueError(error_message)
        self.__num_clusters = num_clusters
        
    def get_num_clusters(self) -> int:
        return self.__num_clusters
        
    def get_coordinates(self) -> typing.Tuple[np.array, np.array]:
        return (self.__x_coord, self.__y_coord)
        
    def get_x_limits(self) -> typing.Tuple[float, float]:
        return (self.__x_min, self.__x_max)

    def get_y_limits(self) -> typing.Tuple[float, float]:
        return (self.__y_min, self.__y_max)
