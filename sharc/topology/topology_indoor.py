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
        self.total_bs_level = self.n_rows*self.n_colums*self.num_cells
        
        self.height = np.empty(0)
        
        
    def calculate_coordinates(self, random_number_gen=np.random.RandomState()):
        """
        Calculates the coordinates of the stations according to the inter-site
        distance parameter. This method is invoked in all snapshots but it can 
        be called only once for the indoor topology. So we set 
        static_base_stations to True to avoid unnecessary calculations.
        """
        if not self.static_base_stations:
            self.static_base_stations = True
            
            x_base = np.array([ (2*k + 1)*self.cell_radius for k in range(self.num_cells)])
            y_base = self.b_d/2*np.ones(self.num_cells)
            
            floor_x = np.empty(0)
            floor_y = np.empty(0)
            for r in range(self.n_rows):
                for c in range(self.n_colums):
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
            for c in range(int(self.n_colums)):
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
    param.n_rows = 2
    param.n_colums = 2
    param.num_floors = 3
    param.street_width = 30
    param.intersite_distance = 30
    param.num_cells = 4
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