# -*- coding: utf-8 -*-
"""
Created on Mon Mar 13 15:37:01 2017

@author: edgar
"""
import numpy as np

from sharc.topology.topology import Topology

class TopologySingleBaseStation(Topology):
    """
    Generates the a single base station centered at (0,0), with azimuth = 0°
    and elevation = -10° wrt x-y plane.
    """
    
    # possible values for base station azimuth and elevation [degrees]
    AZIMUTH = 0
    ELEVATION = -10
    ALLOWED_NUM_CLUSTERS = [1, 2]
    

    def __init__(self, cell_radius: float, num_clusters: int):
        """
        Constructor method that sets the object attributes
        
        Parameters
        ----------
            cell_radius : radius of the cell
        """
        if num_clusters not in TopologySingleBaseStation.ALLOWED_NUM_CLUSTERS:
            error_message = "invalid number of clusters ({})".format(num_clusters)
            raise ValueError(error_message)   
            
        intersite_distance = 2*cell_radius
        super().__init__(intersite_distance, cell_radius)
        self.num_clusters = num_clusters
        
        
    def calculate_coordinates(self):
        """
        Defines the coordinates of the base stations.
        """        
        if not self.static_base_stations:
            self.static_base_stations = True
            if self.num_clusters == 1:
                self.x = np.array([0])
                self.y = np.array([0])
                self.azimuth = TopologySingleBaseStation.AZIMUTH
                self.elevation = TopologySingleBaseStation.ELEVATION
                self.num_base_stations = 1
            elif self.num_clusters == 2:
                self.x = np.array([0, 0])
                self.y = np.array([0, self.intersite_distance])
                self.azimuth = TopologySingleBaseStation.AZIMUTH*np.ones(2)
                self.elevation = TopologySingleBaseStation.ELEVATION*np.ones(2)
                self.num_base_stations = 2             

 