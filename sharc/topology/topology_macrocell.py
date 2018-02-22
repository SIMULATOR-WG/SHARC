# -*- coding: utf-8 -*-
"""
Created on Tue Feb 14 12:51:22 2017

@author: edgar
"""

from sharc.topology.topology import Topology
import matplotlib.pyplot as plt
import matplotlib.axes

import math
import numpy as np

class TopologyMacrocell(Topology):
    """
    Generates the coordinates of the sites based on the macrocell network
    topology.
    """

    # possible values for base station azimuth and elevation [degrees]
    AZIMUTH = [60, 180, 300]
    ELEVATION = -10

    ALLOWED_NUM_CLUSTERS = [1, 7]

    def __init__(self, intersite_distance: float, num_clusters: int):
        """
        Constructor method that sets the parameters and already calls the
        calculation methods.

        Parameters
        ----------
            intersite_distance : Distance between two sites
            num_clusters : Number of clusters, should be 1 or 7
        """
        if num_clusters not in TopologyMacrocell.ALLOWED_NUM_CLUSTERS:
            error_message = "invalid number of clusters ({})".format(num_clusters)
            raise ValueError(error_message)

        cell_radius = intersite_distance*2/3
        super().__init__(intersite_distance, cell_radius)
        self.num_clusters = num_clusters

    def calculate_coordinates(self, random_number_gen=np.random.RandomState()):
        """
        Calculates the coordinates of the stations according to the inter-site
        distance parameter. This method is invoked in all snapshots but it can
        be called only once for the macro cell topology. So we set
        static_base_stations to True to avoid unnecessary calculations.
        """
        if not self.static_base_stations:
            self.static_base_stations = True

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
                    self.y = np.concatenate((self.y, y_central + ys))

            self.x = np.repeat(self.x, 3)
            self.y = np.repeat(self.y, 3)
            self.azimuth = np.tile(self.AZIMUTH, 19*self.num_clusters)
            self.elevation = np.tile(self.ELEVATION, 3*19*self.num_clusters)

            # In the end, we have to update the number of base stations
            self.num_base_stations = len(self.x)
            
            self.indoor = np.zeros(self.num_base_stations, dtype = bool)
                
    def plot(self, ax: matplotlib.axes.Axes):
        # create the hexagons
        r = self.intersite_distance/3
        for x, y, az in zip(self.x, self.y, self.azimuth):
            se = list([[x,y]])
            angle = int(az - 60)
            for a in range(6):
                se.extend([[se[-1][0] + r*math.cos(math.radians(angle)), se[-1][1] + r*math.sin(math.radians(angle))]])
                angle += 60
            sector = plt.Polygon(se, fill=None, edgecolor='k')
            ax.add_patch(sector)

        # macro cell base stations
        ax.scatter(self.x, self.y, color='k', edgecolor="k", linewidth=4, label="Macro cell")


if __name__ == '__main__':
    intersite_distance = 500
    num_clusters = 1
    topology = TopologyMacrocell(intersite_distance, num_clusters)
    topology.calculate_coordinates()

    fig = plt.figure(figsize=(8,8), facecolor='w', edgecolor='k')  # create a figure object
    ax = fig.add_subplot(1, 1, 1)  # create an axes object in the figure

    topology.plot(ax)

    plt.axis('image')
    plt.title("Macro cell topology")
    plt.xlabel("x-coordinate [m]")
    plt.ylabel("y-coordinate [m]")
    plt.tight_layout()
    plt.show()


