"""
Created on Tue Dec 01 11:50:00 2020

@author: Luciano Camilo
"""

import unittest
import numpy as np
from sharc.topology.topology_hibs import TopologyHIBS


class TopologyHIBSTest(unittest.TestCase):

    def setUp(self):
        self.azimuth3 = '60,180,300'
        self.elevation3 = '-90,-90,-90'
        self.azimuth7 = '0,0,60,120,180,240,300'
        self.elevation7 = '-90,-23,-23,-23,-23,-23,-23'
        self.azimuth19 = '0,15,30,45,75,90,105,135,150,165,195,210,225,255,270,285,315,330,345'
        self.elevation19 = '-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30,-30'
        self.num_clusters = 1

    def test_system1(self):

        # Test System 1 3 Sectors
        # Cell Radius 100km

        cell_radius = 100000
        intersite_distance = cell_radius * np.sqrt(3)
        num_sectors = 3
        bs_height = 20000

        topology = TopologyHIBS(intersite_distance, cell_radius, self.num_clusters, num_sectors, bs_height,
                                self.azimuth3, self.azimuth7, self.azimuth19, self.elevation3, self.elevation7,
                                self.elevation19)
        topology.calculate_coordinates()

        # check HIBS Area Radius em Intersite Distance
        self.assertAlmostEqual(topology.intersite_distance, 173205, places=0)
        self.assertEqual(topology.cell_radius, 100000)

        # check HIBS BS Azimuth and Downtilt
        self.assertListEqual(list(topology.azimuth), [60, 180, 300])
        self.assertListEqual(list(topology.elevation), [-90, -90, -90])

        # check HIBS Platform height
        self.assertEqual(topology.height, 20000)
        # check the number of base stations
        self.assertEqual(len(topology.azimuth), topology.num_base_stations)

        # Test System 1 7 Sectors
        # Cell Radius 100km

        cell_radius = 100000
        intersite_distance = cell_radius * np.sqrt(3)
        num_sectors = 7
        bs_height = 20000

        topology = TopologyHIBS(intersite_distance, cell_radius, self.num_clusters, num_sectors, bs_height,
                                self.azimuth3,
                                self.azimuth7,
                                self.azimuth19, self.elevation3, self.elevation7, self.elevation19)
        topology.calculate_coordinates()

        # check HIBS Area Radius em Intersite Distance
        self.assertAlmostEqual(topology.intersite_distance, 173205, places=0)
        self.assertEqual(topology.cell_radius, 100000)

        # check HIBS BS Azimuth and Downtilt
        self.assertListEqual(list(topology.azimuth), [0, 0, 60, 120, 180, 240, 300])
        self.assertListEqual(list(topology.elevation), [-90, -23, -23, -23, -23, -23, -23])

        # check HIBS Platform height
        self.assertEqual(topology.height, 20000)
        # check the number of base stations
        self.assertEqual(len(topology.azimuth), topology.num_base_stations)

        # Test System 1 7 Sectors
        # Cell Radius 90km

        cell_radius = 90000
        intersite_distance = cell_radius * np.sqrt(3)
        num_sectors = 7
        bs_height = 20000

        topology = TopologyHIBS(intersite_distance, cell_radius, self.num_clusters, num_sectors, bs_height,
                                self.azimuth3,
                                self.azimuth7,
                                self.azimuth19, self.elevation3, self.elevation7, self.elevation19)
        topology.calculate_coordinates()

        # check HIBS Area Radius em Intersite Distance
        self.assertAlmostEqual(topology.intersite_distance, 155884.57, places=2)
        self.assertEqual(topology.cell_radius, 90000)

        # check HIBS BS Azimuth and Downtilt
        self.assertListEqual(list(topology.azimuth), [0, 0, 60, 120, 180, 240, 300])
        self.assertListEqual(list(topology.elevation), [-90, -23, -23, -23, -23, -23, -23])

        # check HIBS Platform height
        self.assertEqual(topology.height, 20000)
        # check the number of base stations
        self.assertEqual(len(topology.azimuth), topology.num_base_stations)

    def test_system2(self):

        # Test System 2 19 Sectors
        # Cell Radius 60km
        cell_radius = 60000
        intersite_distance = cell_radius * np.sqrt(3)
        num_sectors = 19
        bs_height = 20000

        topology = TopologyHIBS(intersite_distance, cell_radius, self.num_clusters, num_sectors, bs_height,
                                self.azimuth3,
                                self.azimuth7,
                                self.azimuth19, self.elevation3, self.elevation7, self.elevation19)
        topology.calculate_coordinates()

        # check HIBS Area Radius em Intersite Distance
        self.assertAlmostEqual(topology.intersite_distance, 103923, places=0)
        self.assertEqual(topology.cell_radius, 60000)

        # check HIBS BS Azimuth and Downtilt
        self.assertListEqual(list(topology.azimuth),
                             [0, 15, 30, 45, 75, 90, 105, 135, 150, 165, 195, 210, 225, 255, 270, 285, 315, 330, 345])
        self.assertListEqual(list(topology.elevation),
                             [-30, -30, -30, -30, -30, -30, -30, -30, -30, -30, -30, -30, -30, -30, -30, -30, -30, -30,
                              -30])

        # check HIBS Platform height
        self.assertEqual(topology.height, 20000)
        # check the number of base stations
        self.assertEqual(len(topology.azimuth), topology.num_base_stations)


if __name__ == '__main__':
    unittest.main()
