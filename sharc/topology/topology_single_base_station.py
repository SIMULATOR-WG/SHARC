# -*- coding: utf-8 -*-
"""
Created on Mon Mar 13 15:37:01 2017

@author: edgar
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.axes
import matplotlib.patches as patches

from sharc.topology.topology import Topology

class TopologySingleBaseStation(Topology):
    """
    Generates the a single base station centered at (0,0), with azimuth = 0°
    and elevation = -10° wrt x-y plane.
    """

    # possible values for base station azimuth and elevation [degrees]
    AZIMUTH = [0, 180]
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

    def calculate_coordinates(self, random_number_gen=np.random.RandomState()):
        """
        Defines the coordinates of the base stations.
        """
        if not self.static_base_stations:
            self.static_base_stations = True
            if self.num_clusters == 1:
                self.x = np.array([0])
                self.y = np.array([0])
                self.azimuth = TopologySingleBaseStation.AZIMUTH[0]*np.ones(1)
                self.elevation = TopologySingleBaseStation.ELEVATION*np.ones(1)
                self.num_base_stations = 1
            elif self.num_clusters == 2:
                self.x = np.array([0, self.intersite_distance])
                self.y = np.array([0, 0])
                self.azimuth = np.array(TopologySingleBaseStation.AZIMUTH)
                self.elevation = TopologySingleBaseStation.ELEVATION*np.ones(2)
                self.num_base_stations = 2
            self.indoor = np.zeros(self.num_base_stations, dtype = bool)


    def plot(self, ax: matplotlib.axes.Axes):
        # plot base station
        plt.scatter(self.x, self.y, color='g', edgecolor="w", linewidth=0.5, label="Hotspot")

        # plot base station coverage area
        for x, y, a in zip(self.x, self.y, self.azimuth):
            pa = patches.Wedge( (x, y), self.cell_radius, a-60, a+60, fill=False,
                               edgecolor="green", linestyle='solid' )
            ax.add_patch(pa)


if __name__ == '__main__':
    cell_radius = 100
    num_clusters = 2
    topology = TopologySingleBaseStation(cell_radius, num_clusters)
    topology.calculate_coordinates()

    fig = plt.figure(figsize=(8,8), facecolor='w', edgecolor='k')  # create a figure object
    ax = fig.add_subplot(1, 1, 1)  # create an axes object in the figure

    topology.plot(ax)

    plt.axis('image')
    plt.title("Single base station topology")
    plt.xlabel("x-coordinate [m]")
    plt.ylabel("y-coordinate [m]")
    plt.tight_layout()
    plt.show()

