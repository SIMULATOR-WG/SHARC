# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 15:45:26 2017

@author: edgar
"""

from sharc.topology.topology import Topology
from sharc.parameters.parameters_indoor import ParametersIndoor
import matplotlib.pyplot as plt
import matplotlib.axes

import numpy as np
from itertools import product

class TopologyIndoor(Topology):
    """
    Generates the coordinates of the sites based on the indoor network
    topology. 
    """

    
    def __init__(self, param: ParametersIndoor):
        """
        Constructor method that sets the parameters.
        
        Parameters
        ----------
            param : parameters of the indoor topology
        """

        # These are the building's width, deep and height
        # They do not change
        self.b_w = 120
        self.b_d = 50
        self.b_h = 3
        
        cell_radius = param.intersite_distance/2
        super().__init__(param.intersite_distance, cell_radius)
        
        self.n_rows = param.n_rows
        self.n_colums = param.n_colums
        self.street_width = param.street_width
        self.ue_indoor_percent = param.ue_indoor_percent
        self.building_class = param.building_class
        self.num_cells = param.num_cells
        self.num_floors = param.num_floors
        if param.num_imt_buildings == 'ALL':
            self.all_buildings = True
            self.num_imt_buildings = self.n_rows*self.n_colums
        else:
            self.all_buildings = False
            self.num_imt_buildings = int(param.num_imt_buildings)
        self.imt_buildings = list()
        self.total_bs_level = self.num_imt_buildings*self.num_cells
        
        self.height = np.empty(0)
        
        
    def calculate_coordinates(self, random_number_gen=np.random.RandomState()):
        """
        Calculates the coordinates of the stations according to the inter-site
        distance parameter. This method is invoked in all snapshots but it can 
        be called only once for the indoor topology. So we set 
        static_base_stations to True to avoid unnecessary calculations.
        """
        if not self.static_base_stations:
            self.reset()
            self.static_base_stations = self.all_buildings
            
            x_base = np.array([ (2*k + 1)*self.cell_radius for k in range(self.num_cells)])
            y_base = self.b_d/2*np.ones(self.num_cells)
            
            # Choose random buildings
            all_buildings = list(product(range(self.n_rows),range(self.n_colums)))
            random_number_gen.shuffle(all_buildings)
            self.imt_buildings = all_buildings[:self.num_imt_buildings]
            
            floor_x = np.empty(0)
            floor_y = np.empty(0)
            for build in self.imt_buildings:
                r = build[0]
                c = build[1]
                floor_x = np.concatenate((floor_x, x_base + c*(self.b_w + self.street_width)))
                floor_y = np.concatenate((floor_y, y_base + r*(self.b_d + self.street_width)))

            for f in range(self.num_floors):
                self.x = np.concatenate((self.x, floor_x))
                self.y = np.concatenate((self.y, floor_y))
                self.height = np.concatenate((self.height,
                                             (f+1)*self.b_h*np.ones_like(floor_x)))

            # In the end, we have to update the number of base stations
            self.num_base_stations = len(self.x)        

            self.azimuth = np.zeros(self.num_base_stations)
            self.elevation = -90*np.ones(self.num_base_stations)
            self.indoor = np.ones(self.num_base_stations, dtype = bool)
            
    def reset(self):
        self.x = np.empty(0)
        self.y = np.empty(0)
        self.height = np.empty(0)
        self.azimuth = np.empty(0)
        self.elevation = np.empty(0)
        self.indoor = np.empty(0)
        self.num_base_stations = -1
        self.static_base_stations = False
    
    def plot(self, ax: matplotlib.axes.Axes, top_view = True):
        if top_view:
            self.plot_top_view(ax)
        else:
            self.plot_side_view(ax)
    
    def plot_top_view(self, ax: matplotlib.axes.Axes):
        # create the building
        for b in range(int(self.num_base_stations/self.num_cells)):
            x_b = self.x[self.num_cells*b] - self.cell_radius
            y_b = self.y[self.num_cells*b] - self.b_d/2
            points = list([[x_b, y_b],
                           [x_b + self.b_w, y_b],
                           [x_b + self.b_w, y_b + self.b_d],
                           [x_b, y_b + + self.b_d]])
            sector = plt.Polygon(points, fill=None, edgecolor='k')
            ax.add_patch(sector) 

            for q in range(8):
                points = list()
                x_b = self.x[self.num_cells*b] - self.cell_radius + q*15
                y_b = self.y[self.num_cells*b] + 10
                points.extend([[x_b, y_b],
                               [x_b + 15, y_b],
                               [x_b + 15, y_b + 15],
                               [x_b, y_b + 15]])
                sector = plt.Polygon(points, fill=None, edgecolor='k')
                ax.add_patch(sector)
                
            for q in range(8):
                points = list()
                x_b = self.x[self.num_cells*b] - self.cell_radius + q*15
                y_b = self.y[self.num_cells*b] - self.b_d/2
                points.extend([[x_b, y_b],
                               [x_b + 15, y_b],
                               [x_b + 15, y_b + 15],
                               [x_b, y_b + 15]])
                sector = plt.Polygon(points, fill=None, edgecolor='k')
                ax.add_patch(sector)                
                
        
        # indoor base stations
        ax.scatter(self.x, self.y, color='k', edgecolor="k", linewidth=2, label="Base station")
        
    def plot_side_view(self, ax: matplotlib.axes.Axes):
        
        # Loop on each floor of each column of buildings
        for f in range(int(self.num_floors)):
            for build in self.imt_buildings:
                c = build[1]
                x_b = self.x[f*self.total_bs_level + c*self.num_cells]  - self.cell_radius 
                z_b = self.height[f*self.total_bs_level + c*self.num_cells]
                points = list([[x_b, z_b],
                               [x_b + self.b_w, z_b],
                               [x_b + self.b_w, z_b - self.b_h],
                               [x_b, z_b - self.b_h]])
                sector = plt.Polygon(points, fill=None, edgecolor='k')
                ax.add_patch(sector) 
        
        ax.scatter(self.x, self.height-0.05, color='k', edgecolor="k", linewidth=2, label="Base station")

if __name__ == '__main__':
    param = ParametersIndoor()
    param.intersite_distance = 20
    param.n_rows = 5
    param.n_colums = 5
    param.num_imt_buildings = 5
    param.num_floors = 3
    param.street_width = 30
    param.intersite_distance = 20
    param.num_cells = 6
    param.ue_indoor_percent = 0.95
    param.building_class = "TRADITIONAL"
    topology = TopologyIndoor(param)
    topology.calculate_coordinates()
    
    # Plot top view
    fig = plt.figure(figsize=(10,8), facecolor='w', edgecolor='k')  # create a figure object
    ax = fig.add_subplot(1, 1, 1)  # create an axes object in the figure
    
    topology.plot(ax)
    
    plt.axis('image') 
    plt.title("Indoor topology")
    plt.xlabel("x-coordinate [m]")
    plt.ylabel("y-coordinate [m]")
    plt.tight_layout()    
    
    axes = plt.gca()
#    axes.set_xlim([-param.street_width, param.n_colums*3*param.intersite_distance + (param.n_colums-1)*param.street_width + param.street_width])
#    axes.set_ylim([-param.street_width, param.n_rows*topology.b_d + (param.n_rows-1)*param.street_width + param.street_width])
    
    plt.show()    
    
    # Plot side view
    fig = plt.figure(figsize=(10,8), facecolor='w', edgecolor='k')  # create a figure object
    ax = fig.add_subplot(1, 1, 1)  # create an axes object in the figure
    
    topology.plot(ax,top_view=False)
    
    plt.title("Indoor topology")
    plt.xlabel("x-coordinate [m]")
    plt.ylabel("z-coordinate [m]")  
    plt.tight_layout()
    
    axes = plt.gca()
    axes.set_ylim((0,3*param.num_floors + 1))
    plt.show()  