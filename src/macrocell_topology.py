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
    
    def __init__(self, intersite_distance: float, num_clusters: int):
        self.__ALLOWED_NUM_CLUSTERS = [1,7]
        self.set_intersite_distance(intersite_distance)
        self.set_num_clusters(num_clusters)
    
    def set_intersite_distance(self, intersite_distance: float):
        self.__intersite_distance = intersite_distance
        
    def get_intersite_distance(self) -> float:
        return self.__intersite_distance       
        
    def set_num_clusters(self, num_clusters: int):
        if num_clusters not in self.__ALLOWED_NUM_CLUSTERS:
            error_message = "invalid number of clusters ({})".format(num_clusters)
            raise ValueError(error_message)
        self.__num_clusters = num_clusters
        
    def get_num_clusters(self) -> int:
        return self.__num_clusters
        
    def get_coordinates(self) -> typing.Tuple[np.array, np.array]:
        d = self.get_intersite_distance()
        h = (d/3)*math.sqrt(3)/2
        # these are the coordinates of the central cluster
        x_central = np.array([0, d, d/2, -d/2, -d, -d/2, 
                         d/2, 2*d, 3*d/2, d, 0, -d, 
                         -3*d/2, -2*d, -3*d/2, -d, 0, d, 3*d/2])
        y_central = np.array([0, 0, 3*h, 3*h, 0, -3*h, 
                         -3*h, 0, 3*h, 6*h, 6*h, 6*h, 
                         3*h, 0, -3*h, -6*h, -6*h, -6*h, -3*h])
        x_coord = np.copy(x_central)
        y_coord = np.copy(y_central)
        # other clusters are calculated by shifting the central cluster
        if self.get_num_clusters() == 7:
            x_shift = np.array([7*d/2, -d/2, -4*d, -7*d/2, d/2, 4*d])
            y_shift = np.array([9*h, 15*h, 6*h, -9*h, -15*h, -6*h])
            for x, y in zip(x_shift, y_shift):
                x_coord = np.concatenate((x_coord, x_central+x))
                y_coord = np.concatenate((y_coord, y_central+y))
        
        return (x_coord, y_coord)