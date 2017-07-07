# -*- coding: utf-8 -*-
"""
Created on Tue May 16 16:59:40 2017

@author: edgar
"""

import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.patches as patches

from shapely.geometry import Polygon

from sharc.topology.topology import Topology
from sharc.topology.topology_macrocell import TopologyMacrocell
from sharc.parameters.parameters_hotspot import ParametersHotspot


class TopologyHotspot(Topology):
    """
    Generates the coordinates of the stations based on the heterogeneous 
    network topology (macro cell with hotspots).
    """
    
    # Possible values for base station azimuth [degrees]. 
    # The value is randomly chosen from this array
    AZIMUTH = [0, 90, 180, 270]   

    # Posible values for base station elevation [degrees]
    ELEVATION = -10 
    
    def __init__(self, param: ParametersHotspot, intersite_distance: float, num_clusters: int):
        """
        Constructor method that sets the parameters and already calls the 
        calculation methods.
        
        Parameters
        ----------
            param : Hotspots parameters
            intersite_distance : Distance between macro cell base stations
            num_clusters : Number of macro cell cluters, should be 1 or 7
        """
        self.param = param
        self.macrocell = TopologyMacrocell(intersite_distance, num_clusters)
        self.macrocell.calculate_coordinates()

        cell_radius = self.param.max_dist_hotspot_ue
        super().__init__(intersite_distance, cell_radius)
        
        
    def calculate_coordinates(self):
        """
        Calculates coordinates of hotspots
        """
        i = 0
        x = np.empty(0)
        y = np.empty(0)
        azimuth = np.empty(0)
        for cell_x, cell_y, cell_azimuth in zip(self.macrocell.x, self.macrocell.y, self.macrocell.azimuth):
            #print("base station #{}".format(i))
            i += 1
            # find the center coordinates of the sector (hexagon)
            macro_cell_x = cell_x + self.macrocell.intersite_distance/3*math.cos(math.radians(cell_azimuth))
            macro_cell_y = cell_y + self.macrocell.intersite_distance/3*math.sin(math.radians(cell_azimuth))
            # generate hotspots center coordinates
            hotspots_validated = False
            while(not hotspots_validated):
                # Hotspots are generated inside an inscribed circle of a regular hexagon (sector).
                # The backoff factor (1.0) controls the overlapping rate between hotspots
                # coverage areas (overlapping of hotspots in different macro cells)
                r = (self.macrocell.intersite_distance/3)*np.sqrt(3)/2 - self.param.max_dist_hotspot_ue/1.0
                hotspot_radius = r*np.random.random(self.param.num_hotspots_per_cell)
                hotspot_angle = 2*np.pi*np.random.random(self.param.num_hotspots_per_cell)
                hotspot_x = hotspot_radius*np.cos(hotspot_angle) + macro_cell_x
                hotspot_y = hotspot_radius*np.sin(hotspot_angle) + macro_cell_y
                hotspot_azimuth = np.random.choice(self.AZIMUTH, self.param.num_hotspots_per_cell)
                # Hotspots within a cell are validated if they do not overlap
                # and if they have the minimum separation distance from macro BS
                hotspots_validated = (not self.overlapping_hotspots(hotspot_x,
                                                                    hotspot_y,
                                                                    hotspot_azimuth,
                                                                    self.cell_radius*np.ones(self.param.num_hotspots_per_cell))) and \
                                          self.validade_min_dist_bs_hotspot(hotspot_x, 
                                                                            hotspot_y, 
                                                                            self.macrocell.x, 
                                                                            self.macrocell.y, 
                                                                            self.param.min_dist_bs_hotspot)
            x = np.concatenate([x, hotspot_x])
            y = np.concatenate([y, hotspot_y])   
            azimuth = np.concatenate([azimuth, hotspot_azimuth])   

        self.x = x
        self.y = y
        self.azimuth = azimuth
        self.elevation = self.ELEVATION*np.ones(len(self.x))
        # In the end, we have to update the number of base stations
        self.num_base_stations = len(self.x)
                
                
    def overlapping_hotspots(self, 
                             x: np.array, 
                             y: np.array,
                             azimuth: np.array,
                             radius: np.array) -> bool:
        """
        Evaluates the spatial relationships among hotspots and checks whether
        hotspots coverage areas intersect.
        
        Parameters
        ----------
            x: x-coordinates of the hotspots
            y: y-coordinates of the hotspots
            azimuth: horizontal angle of the hotspots (orientation)
            radius: radius of the coverage area of the hotspots
            
        Returns
        -------
            True if there is intersection between any two hotspots
        """
        # Each hotspot coverage area corresponds to a Polygon object
        polygons = list()
        for x, y, azimuth, r in zip(x, y, azimuth, radius):
            points = list()
            points.append((x, y))
            azimuth_values = np.linspace(-60, 60, 25)
            for a in range(len(azimuth_values)):
                points.append((x + r*math.cos(np.radians(azimuth + azimuth_values[a])), y + r*math.sin(np.radians(azimuth + azimuth_values[a]))))
            polygons.append(Polygon(points))
        
        # Check if there is overlapping between any of the hotspots coverage
        # areas. In other words, check if any polygons intersect 
        for p in range(len(polygons)-1):
            for pi in range(p+1, len(polygons)):
                overlapping = polygons[p].intersects(polygons[pi])
                if overlapping:
                    # An intersection was found! We stop here because we do not
                    # need to check other combinations
                    return True
        
        # If this point is reached, then there is no intersection between polygons
        return False
        
        
    def validade_min_dist_hotspots(self, 
                                   hotspot_x: np.array, 
                                   hotspot_y: np.array, 
                                   min_dist_hotspots: float) -> bool:
        """
        Checks minimum 2D distance between two hotspots. Currently not used.
        
        Returns
        -------
        out : bool
            True if hotspots coordinates meets the minimum 2D distance between 
            any two hotspots
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
        hotspots.
        
        Returns
        -------
        out : bool
            True if hotspots coordinates meets the minimum 2D distance between 
            macro cell base stations and hotspots
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



def plot_topology(topology: TopologyHotspot):
    fig = plt.figure(figsize=(8,8), facecolor='w', edgecolor='k')
    ax = fig.gca()
    
    r = topology.macrocell.intersite_distance/3
    for x, y, az in zip(topology.macrocell.x, topology.macrocell.y, topology.macrocell.azimuth):
        se = list([[x,y]])
        angle = int(az - 60)
        for a in range(6):
            se.extend([[se[-1][0] + r*math.cos(math.radians(angle)), se[-1][1] + r*math.sin(math.radians(angle))]])
            angle += 60
        sector = plt.Polygon(se, fill=None, edgecolor='k')
        ax.add_patch(sector)
    
    # plot small cells
    plt.scatter(topology.x, topology.y, color='g', edgecolor="w", linewidth=0.5, label="Hotspot")
    
    # plot small cells coverage area
    for x, y, a in zip(topology.x, topology.y, topology.azimuth):
        pa = patches.Wedge( (x, y), topology.cell_radius, a-60, a+60, fill=False, 
                           edgecolor="green", linestyle='solid' )
        ax.add_patch(pa)

    
    # macro cell base stations
    plt.scatter(topology.macrocell.x, topology.macrocell.y, color='k', edgecolor="k", linewidth=4, label="Macro cell")

    # sector centers
    #plt.scatter(-sector_y, sector_x, color='g', edgecolor="g")
    
    # plot macro cell coverage area
#    r = (topology.macrocell.intersite_distance/3)*math.sqrt(3)/2 - topology.param.max_dist_hotspot_ue/2
#    for x, y, az in zip(topology.macrocell.x, topology.macrocell.y, topology.macrocell.azimuth):
#        # find the center coordinates of the sector (hexagon)
#        mx = x + topology.macrocell.intersite_distance/3*math.cos(math.radians(az))
#        my = y + topology.macrocell.intersite_distance/3*math.sin(math.radians(az))
#        circ = plt.Circle((mx, my), radius=r, color='b', fill=False, linewidth=0.5)
#        ax.add_patch(circ)    
    
    plt.axis('image') 
    plt.title("Macro cell topology with hotspots")
    plt.xlabel("x-coordinate [m]")
    plt.ylabel("y-coordinate [m]")
    #plt.xlim((-3000, 3000))
    #plt.ylim((-3000, 3000))                
    plt.legend(loc="upper left", scatterpoints=1)
    plt.tight_layout()    
    plt.show()


if __name__ == '__main__':
    param = ParametersHotspot()
    param.num_hotspots_per_cell = 2
    param.max_dist_hotspot_ue = 60
    param.min_dist_hotspot_ue = 5
    param.min_dist_bs_hotspot = 100
    param.min_dist_hotspots = 2*param.max_dist_hotspot_ue    

    intersite_distance = 500
    num_clusters = 1
    topology = TopologyHotspot(param, intersite_distance, num_clusters)
    topology.calculate_coordinates()
    plot_topology(topology)

