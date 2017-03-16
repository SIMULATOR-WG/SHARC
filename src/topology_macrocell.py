# -*- coding: utf-8 -*-
"""
Created on Tue Feb 14 12:51:22 2017

@author: edgar
"""

from topology import Topology

import math
import numpy as np

class TopologyMacrocell(Topology):
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
        cell_radius = intersite_distance*2/3
        super(TopologyMacrocell, self).__init__(intersite_distance, 
                                                cell_radius, num_clusters)

    def _calculate_coordinates(self):
        """
        Calculates the coordinates of the stations according to the inter-site
        distance parameter.
        """
        d = self.intersite_distance
        h = (d/3)*math.sqrt(3)/2
        # these are the coordinates of the central cluster
        x_central = np.array([0, d, d/2, -d/2, -d, -d/2, 
                         d/2, 2*d, 3*d/2, d, 0, -d, 
                         -3*d/2, -2*d, -3*d/2, -d, 0, d, 3*d/2])
        y_central = np.array([0, 0, 3*h, 3*h, 0, -3*h, 
                         -3*h, 0, 3*h, 6*h, 6*h, 6*h, 
                         3*h, 0, -3*h, -6*h, -6*h, -6*h, -3*h])
        self.x = np.copy(x_central)
        self.y = np.copy(y_central)
        # other clusters are calculated by shifting the central cluster
        if self.num_clusters == 7:
            x_shift = np.array([7*d/2, -d/2, -4*d, -7*d/2, d/2, 4*d])
            y_shift = np.array([9*h, 15*h, 6*h, -9*h, -15*h, -6*h])
            for xs, ys in zip(x_shift, y_shift):
                self.x = np.concatenate((self.x, x_central + xs))
                self.y = np.concatenate((self.x, y_central + ys))    
    
    def _calculate_limits(self):
        """
        Calculates the coordinates of the scenario's borders
        """
        self.x_min = np.min(self.x) - self.cell_radius
        self.x_max = np.max(self.x) + self.cell_radius
        self.y_min = np.min(self.y) - self.cell_radius
        self.y_max = np.max(self.y) + self.cell_radius

    @Topology.num_clusters.setter
    def num_clusters(self, value):
        """
        Override the definition and check if number of clusters is valid
        """
        if value not in self.__ALLOWED_NUM_CLUSTERS:
            error_message = "invalid number of clusters ({})".format(value)
            raise ValueError(error_message)
        Topology.num_clusters.fset(self, value)
        
    @Topology.cell_radius.setter
    def cell_radius(self, value):
        """
        When cell radius changes, intersite distance also has to change
        """
        self._intersite_distance = value*3/2
        Topology.cell_radius.fset(self, value)
        
    @Topology.intersite_distance.setter
    def intersite_distance(self, value):
        """
        When intersite distance changes, cell radius also has to change
        """
        self._cell_radius = value*2/3
        Topology.intersite_distance.fset(self, value)
      
