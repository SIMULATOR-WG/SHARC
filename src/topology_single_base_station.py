# -*- coding: utf-8 -*-
"""
Created on Mon Mar 13 15:37:01 2017

@author: edgar
"""

from topology import Topology

class TopologySingleBaseStation(Topology):
    """
    Generates the a single base station centered at (0,0).
    """
    
    def __init__(self, cell_radius: float):
        """
        Constructor method that sets the parameters and already calls the 
        calculation methods.
        
        Parameters
        ----------
            cell_radius : radius of the cell
        """
        # Actually, intersite distance is not used
        intersite_distance = 2*cell_radius
        # Naturally, number of clusters has to be equal to 1
        num_clusters = 1
        super(TopologySingleBaseStation, self).__init__(intersite_distance, 
                                                cell_radius, num_clusters)
        
    def _calculate_coordinates(self):
        """
        Defines the coordinates of the station.
        """        
        self.x = 0
        self.y = 0
        
    def _calculate_limits(self):
        """
        Calculates the coordinates of the scenario's borders
        """
        self.x_min = self.x - self.cell_radius
        self.x_max = self.x + self.cell_radius
        self.y_min = self.y - self.cell_radius
        self.y_max = self.y + self.cell_radius

    @Topology.num_clusters.setter
    def num_clusters(self, value):
        """
        The single base station topology supports only one cluster. Any attempt
        to set this value will raise an error.
        """
        error_message = "Cannot set number of clusters in single base station network topology"
        raise NotImplementedError(error_message)
        
    @Topology.intersite_distance.setter
    def intersite_distance(self, value):
        """
        The intersite distance attribute is meaningless in the single base 
        station topology. Any attempt to set this value will raise an error.
        """
        error_message = "Cannot set intersite distance in single base station network topology"
        raise NotImplementedError(error_message)
 