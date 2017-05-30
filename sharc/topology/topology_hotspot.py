# -*- coding: utf-8 -*-
"""
Created on Tue May 16 16:59:40 2017

@author: edgar
"""

import numpy as np
import math

#from sharc.topology.topology import Topology
from sharc.topology.topology_macrocell import TopologyMacrocell
from sharc.parameters.parameters_hotspot import ParametersHotspot


class TopologyHotspot(object):
    """
    Generates the coordinates of the stations based on the heterogeneous 
    network topology (macro cell with hotspots).
    """
    
    def __init__(self, param: ParametersHotspot, intersite_distance: float, num_clusters: int):
        """
        Constructor method that sets the parameters and already calls the 
        calculation methods.
        
        Parameters
        ----------
            intersite_distance : Distance between macro cell base stations
            num_clusters : Number of macro cell cluters, should be 1 or 7
        """

        self.topology_macrocell = TopologyMacrocell(intersite_distance, num_clusters)
        self.param = param
        
        self.x = np.empty(0)
        self.y = np.empty(0)    
        
        self.hotspot_x = np.empty(0)
        self.hotspot_y = np.empty(0)
        
        
    def calculate_coordinates(self):
        """
        Calculates coordinates of macrocell and hotspots base stations
        """
        i = 0
        for cell_x, cell_y in zip(self.topology_macrocell.x, self.topology_macrocell.y):
            azimuth = np.radians([60, 180, 300])
            print("site #{}".format(i))
            i += 1
            for a in azimuth:
                # find the center coordinates of the sector (hexagon)
                macro_cell_x = cell_x + self.topology_macrocell.intersite_distance/3*math.cos(a)
                macro_cell_y = cell_y + self.topology_macrocell.intersite_distance/3*math.sin(a)
                #macro_cell_x = cell_x
                #macro_cell_y = cell_y
                # x = (x_max - x_min)*np.random.random(num_samples) + x_min
                # generate hotspots center coordinates
                hotspots_validated = False
                while(not hotspots_validated):
                    # hotspots are generated inside a inscribed circle of a regular hexagon (sector)
                    r = (self.topology_macrocell.intersite_distance/3)*np.sqrt(3)/2 - self.param.max_dist_ue_hotspot
                    hotspot_radius = r*np.random.random(self.param.num_hotspots_per_cell)
                    hotspot_angle = 2*np.pi*np.random.random(self.param.num_hotspots_per_cell)
                    hotspot_x = hotspot_radius*np.cos(hotspot_angle) + macro_cell_x
                    hotspot_y = hotspot_radius*np.sin(hotspot_angle) + macro_cell_y
                    hotspots_validated = \
                        self.validade_min_dist_hotspots(hotspot_x, 
                                                        hotspot_y, 
                                                        self.param.min_dist_hotspots) and \
                        self.validade_min_dist_bs_hotspot(hotspot_x, 
                                                          hotspot_y, 
                                                          self.topology_macrocell.x, 
                                                          self.topology_macrocell.y, 
                                                          self.param.min_dist_bs_hotspot)
                self.hotspot_x = np.concatenate([self.hotspot_x, hotspot_x])
                self.hotspot_y = np.concatenate([self.hotspot_y, hotspot_y])                    
                # generate base station coordinates within hotspot area
                for x, y in zip(hotspot_x, hotspot_y):
                    stations_validated = False
                    while(not stations_validated):
                        station_radius = self.param.max_dist_station_hotspot*np.random.random(self.param.num_stations_per_hotspots)
                        station_angle = 2*np.pi*np.random.random(self.param.num_stations_per_hotspots)
                        station_x = station_radius*np.cos(station_angle) + x
                        station_y = station_radius*np.sin(station_angle) + y
                        stations_validated = \
                            self.validade_min_dist_hotspots(station_x, 
                                                            station_y, 
                                                            self.param.min_dist_stations)
                    self.x = np.concatenate([self.x, station_x])
                    self.y = np.concatenate([self.y, station_y])
                
                
    def validade_min_dist_hotspots(self, 
                                   hotspot_x: np.array, 
                                   hotspot_y: np.array, 
                                   min_dist_hotspots: float) -> bool:
        """
        Checks minimum 2D distance between two hotspot centers.
        
        Returns
        -------
        out : bool
            True if hotspots coordinates meets the minimum 2D distance between 
            any two hotspot centers
        """
        # Here we have a 2D matrix whose values indicates the distance between
        # the hotspots. The diagonal elements are obviously equal to zero
        distance = np.sqrt((hotspot_x - hotspot_x.reshape((-1, 1)))**2 + 
                           (hotspot_y - hotspot_y.reshape((-1, 1)))**2)
        num_hotpots = len(hotspot_x)
        # count the number of values that are less than the minimum distance and
        # return true if this value is equal os less than the number of hotspots.
        # In other words, it returns True if only diagonal elements are less
        # than the minimum distance
        occ = np.where(distance < min_dist_hotspots)[0]
        return len(occ) == num_hotpots


    def validade_min_dist_bs_hotspot(self, 
                                     hotspot_x: np.array, 
                                     hotspot_y: np.array, 
                                     macrocell_x: np.array, 
                                     macrocell_y: np.array, 
                                     min_dist_bs_hotspot: float) -> bool:
        """
        Checks minimum 2D distance between macro cell base stations and 
        hotspot centers.
        
        Returns
        -------
        out : bool
            True if hotspots coordinates meets the minimum 2D distance between 
            macro cell base stations and hotspot centers
        """
        # Here we have a 2D matrix whose values indicates the distance between
        # base station and hotspots. In this matrix, each line corresponds to 
        # a macro cell base station and each column corresponds to a hotspot
        distance = np.sqrt((hotspot_x - macrocell_x.reshape((-1, 1)))**2 + 
                           (hotspot_y - macrocell_y.reshape((-1, 1)))**2)
        # count the number of values that are less than the minimum distance and
        # return true if any value is equal os less than minimum 2D distance 
        # between macro cell base stations and hotspot centers
        occ = np.where(distance < min_dist_bs_hotspot)[0]
        return len(occ) == 0


    def validade_min_dist_stations(self, 
                                   station_x: np.array, 
                                   station_y: np.array,  
                                   min_dist_stations: float) -> bool:
        """
        Checks minimum 2D distance between two stations (small cells) in a 
        hotspot.
        
        Returns
        -------
        out : bool
            True if stations' coordinates meets the minimum 2D distance separation 
        """
        # Here we have a 2D matrix whose values indicates the distance between
        # the stations. The diagonal elements are obviously equal to zero
        distance = np.sqrt((station_x - station_x.reshape((-1, 1)))**2 + 
                           (station_y - station_y.reshape((-1, 1)))**2)
        num_stations = len(station_x)
        # count the number of values that are less than the minimum distance and
        # return true if this value is equal os less than the number of stations.
        # In other words, it returns True if only diagonal elements are less
        # than the minimum distance
        occ = np.where(distance < min_dist_stations)[0]
        return len(occ) == num_stations