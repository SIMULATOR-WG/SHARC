# -*- coding: utf-8 -*-
"""
Created on Mon Mar 13 15:37:01 2017

@author: edgar
"""
import numpy as np

from sharc.topology.topology import Topology

class TopologySingleBaseStation(Topology):
    """
    Generates the a single base station centered at (0,0).
    """
    
    def __init__(self, cell_radius: float, num_clusters: int):
        """
        Constructor method that sets the parameters and already calls the 
        calculation methods.
        
        Parameters
        ----------
            cell_radius : radius of the cell
        """
        allowed_num_clusters = [1,2]
        # Actually, intersite distance is not used
        intersite_distance = 2*cell_radius
        super().__init__(intersite_distance, cell_radius, num_clusters, allowed_num_clusters)
        
    def calculate_coordinates(self):
        """
        Defines the coordinates of the station.
        """        
        if self.num_clusters == 1:
            self.x = np.array([0])
            self.y = np.array([0])
        elif self.num_clusters == 2:
            self.x = np.array([-self.cell_radius, self.cell_radius])
            self.y = np.array([0, 0])
        else:
            error_message = "invalid number of clusters ({})".format(self.num_clusters)
            raise ValueError(error_message)             
        

 