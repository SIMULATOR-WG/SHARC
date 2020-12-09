# -*- coding: utf-8 -*-
"""
Created on Thu Nov 12 21:26:22 2020

@author: Luciano Camilo
"""

from sharc.topology.topology import Topology
import matplotlib.pyplot as plt
import matplotlib.axes
import math
import numpy as np


class TopologyHIBS(Topology):
    """
    ____________________________________________________________________________________
    Generates the coordinates of the HIBS sites based on the macrocell network topology.
    ____________________________________________________________________________________
    """

    ALLOWED_NUM_CLUSTERS = [1]

    def __init__(self,
                 intersite_distance_hibs: float,
                 cell_radius_hibs: int,
                 num_clusters_hibs: 1,
                 number_of_sectors=None,
                 bs_height_hibs=20000,
                 azimuth3_hibs=None,
                 azimuth7_hibs=None,
                 azimuth19_hibs=None,
                 elevation3_hibs=None,
                 elevation7_hibs=None,
                 elevation19_hibs=None):
        """
        Constructor method that sets the parameters and already calls the
        calculation methods.

        Parameters
        ----------
            intersite_distance_hibs : Distance between two sites
            num_clusters_hibs : Number of clusters, should be 1 or 7
        """

        if num_clusters_hibs not in TopologyHIBS.ALLOWED_NUM_CLUSTERS:
            error_message = "invalid number of clusters ({})".format(num_clusters_hibs)
            raise ValueError(error_message)

        if number_of_sectors not in [1, 3, 7, 19]:
            error_message = "invalid number of sectors ({})".format(num_clusters_hibs)
            raise ValueError(error_message)

        # cell_radius = intersite_distance / np.sqrt(3)

        super().__init__(intersite_distance_hibs, cell_radius_hibs)
        self.cell_radius = cell_radius_hibs
        self.intersite_distance = cell_radius_hibs * np.sqrt(3)
        self.num_clusters = num_clusters_hibs
        self.num_sectors = number_of_sectors
        self.height = bs_height_hibs
        self.azimuth3 = azimuth3_hibs
        self.azimuth7 = azimuth7_hibs
        self.azimuth19 = azimuth19_hibs
        self.elevation3 = elevation3_hibs
        self.elevation7 = elevation7_hibs
        self.elevation19 = elevation19_hibs

    def calculate_coordinates(self, random_number_gen=np.random.RandomState()):
        """
        Calculates the coordinates of the stations according to the inter-site
        distance parameter. This method is invoked in all snapshots but it can
        be called only once for the macro cell topology. So we set
        static_base_stations to True to avoid unnecessary calculations.
        """

        global azimuth, elevation
        if not self.static_base_stations:

            self.static_base_stations = True

            aux_azi3 = self.azimuth3.split(',')
            self.azimuth3 = [int(i) for i in aux_azi3]

            aux_azi7 = self.azimuth7.split(',')
            self.azimuth7 = [int(i) for i in aux_azi7]

            aux_azi19 = self.azimuth19.split(',')
            self.azimuth19 = [int(i) for i in aux_azi19]

            aux_ele3 = self.elevation3.split(',')
            self.elevation3 = [int(i) for i in aux_ele3]

            aux_ele7 = self.elevation7.split(',')
            self.elevation7 = [int(i) for i in aux_ele7]

            aux_ele19 = self.elevation19.split(',')
            self.elevation19 = [int(i) for i in aux_ele19]

            # d = self.intersite_distance / 2
            # h = self.cell_radius

            # these are the coordinates of the central cluster
            x_central = np.array([0])
            y_central = np.array([0])

            self.x = np.copy(x_central)
            self.y = np.copy(y_central)
            """
            ___________________________________________________________________________________________________________
            Fix the azimuth and elevation values for each cell configuration according ITU document WP 5D 237-E.
            ___________________________________________________________________________________________________________
            """
            if self.num_sectors == 1:  # 1 Sector (Test System)
                # AZIMUTH = [60, 180, 300]
                azimuth = [0]
                # ELEVATION = [-90, -90, -90]
                elevation = [-90]

            if self.num_sectors == 3:  # 3 Sectors (Possible System 1)
                # AZIMUTH = [60, 180, 300]
                azimuth = self.azimuth3
                # ELEVATION = [-90, -90, -90]
                elevation = self.elevation3

            elif self.num_sectors == 7:  # 7 Sectors  (System 1)
                # AZIMUTH = [0, 30, 60, 90, 120, 180, 240]
                azimuth = self.azimuth7
                # ELEVATION = [-90, -23, -23, -23, -23, -23, -23]
                elevation = self.elevation7

            elif self.num_sectors == 19:   # 19 Sectors (System 2)
                # AZIMUTH = [0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 165, 180, 195, 225, 240, 255, 285, 315, 345]
                azimuth = self.azimuth19
                # ELEVATION = [-30, -30, -30, -30, -30, -30, -30, -30, -30, -30, -30, -30, -30, -30, -30, -30, -30,
                #             -30, -30 ]
                elevation = self.elevation19

            self.x = np.repeat(self.x, self.num_sectors)
            self.y = np.repeat(self.y, self.num_sectors)

            if self.num_sectors in [1, 3, 7, 19]:
                self.azimuth = np.array(azimuth)
                self.elevation = np.array(elevation)

            # In the end, we have to update the number of base stations
            self.num_base_stations = len(self.x)

            self.indoor = np.zeros(self.num_base_stations, dtype=bool)

    def plot(self, axis: matplotlib.axes.Axes):

        # create the hexagons for 1 Sectors
        global azimuth, r, angle
        if self.num_sectors == 1:
            r = self.cell_radius
            azimuth = np.zeros(self.x.size)

        # create the hexagons for 3 Sectors
        elif self.num_sectors == 3:
            r = self.intersite_distance / np.sqrt(3)
            azimuth = self.azimuth

        # create the hexagon for 7 Sectors
        elif self.num_sectors == 7:
            r = self.cell_radius
            azimuth = np.zeros(self.x.size)

        # create the dodecagon for 19 Sectors
        elif self.num_sectors == 19:
            r = self.intersite_distance / np.sqrt(3)
            azimuth = np.zeros(self.x.size)

        # create the hexagons for 3 Sectors
        for x, y, az in zip(self.x, self.y, azimuth):

            if self.num_sectors == 1:
                x = x - self.intersite_distance / 2
                y = y - r / 2
                angle = int(az - 30)

            elif self.num_sectors == 3:
                angle = int(az - 30)

        # create the hexagon for 7 Sectors
            elif self.num_sectors == 7:
                x = x - self.intersite_distance / np.sqrt(3)
                # y = y - r /2
                y = y
                angle = int(az - 60)

        # create the dodecagon for 19 Sectors
            elif self.num_sectors == 19:
                x = x
                y = y
                angle = int(az)

            se = list([[x, y]])

            # plot polygon - 1 Sectors
            if self.num_sectors == 1:
                for a in range(7):
                    se.extend(
                        [[se[-1][0] + r * math.cos(math.radians(angle)), se[-1][1] +
                          r * math.sin(math.radians(angle))]])
                    angle += 60
                sector = plt.Polygon(se, fill=None, edgecolor='k')
                axis.add_patch(sector)

                # plot polygon - 3 Sectors
            elif self.num_sectors == 3:
                for a in range(6):
                    se.extend(
                        [[se[-1][0] + r / 2 * math.cos(math.radians(angle)), se[-1][1] +
                          r / 2 * math.sin(math.radians(angle))]])
                    angle += 60
                sector = plt.Polygon(se, fill=None, edgecolor='blue', linewidth=1, alpha=1)
                axis.add_patch(sector)

            # plot polygon - 7 Sectors
            elif self.num_sectors == 7:
                for a in range(7):
                    se.extend(
                        [[se[-1][0] + r * math.cos(math.radians(angle)), se[-1][1] +
                          r * math.sin(math.radians(angle))]])
                    angle += 60
                sector = plt.Polygon(se, fill=None, edgecolor='k')
                axis.add_patch(sector)

                # plot polygon - 19 Sectors
            elif self.num_sectors == 19:
                for a in range(25):
                    se.extend(
                        [[se[0][0] + r * math.cos(math.radians(angle)),
                          se[0][0] + r * math.sin(math.radians(angle))]])
                    angle += 15
                sector = plt.Polygon(se[::-2], fill=None, edgecolor='k')
                axis.add_patch(sector)

        # macro cell base stations
        axis.scatter(self.x, self.y, s=200, marker='v', c='k', edgecolor='k', linewidth=1, alpha=1,
                     label="HIBs platforms")


if __name__ == '__main__':

    cell_radius = 100000
    intersite_distance = cell_radius * np.sqrt(3)
    num_clusters = 1
    num_sectors = 1
    bs_height = 20000
    azimuth3 = '60,180,300'
    elevation3 = '-90,-90,-90'
    azimuth7 = '0,0,60,120,180,240,300'
    elevation7 = '-90,-23,-23,-23,-23,-23,-23'
    azimuth19 = '0,15,30,45,75,90,105,135,150,165,195,210,225,255,270,285,315,330,345'
    elevation19 = '-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30'

    topology = TopologyHIBS(intersite_distance, cell_radius, num_clusters, num_sectors, bs_height, azimuth3, azimuth7,
                            azimuth19, elevation3, elevation7, elevation19)
    topology.calculate_coordinates()

    fig = plt.figure(figsize=(8, 8), facecolor='w', edgecolor='k')  # create a figure object
    ax = fig.add_subplot(1, 1, 1)  # create an axes object in the figure

    topology.plot(ax)

    plt.axis('image')
    plt.title("HIBS System 1 Topology - 7 Sector")

    plt.xlabel("x-coordinate [km]")
    plt.ylabel("y-coordinate [km]")
    plt.legend()
    plt.tight_layout()
    plt.show()
