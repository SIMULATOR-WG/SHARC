# -*- coding: utf-8 -*-
"""
Created on Mon Mar 13 17:21:44 2017

@author: edgar
"""

import unittest
import numpy as np
import numpy.testing as npt

from sharc.topology.topology_single_base_station import TopologySingleBaseStation

class TopologySingleBaseStationTest(unittest.TestCase):
    
    def setUp(self):
        pass
        
    def test_coordinates(self):
        cell_radius = 1500
        num_clusters = 1
        topology = TopologySingleBaseStation(cell_radius, num_clusters)
        topology.calculate_coordinates()
        npt.assert_equal(topology.x, np.array([0]))
        npt.assert_equal(topology.y, np.array([0]))

        cell_radius = 1500
        num_clusters = 2
        topology = TopologySingleBaseStation(cell_radius, num_clusters)
        topology.calculate_coordinates()
        npt.assert_equal(topology.x, np.array([0, 2*cell_radius]))
        npt.assert_equal(topology.y, np.array([0, 0]))

        
    def test_num_clusters(self):
        cell_radius = 1500
        num_clusters = 1
        topology = TopologySingleBaseStation(cell_radius, num_clusters)
        topology.calculate_coordinates()
        self.assertEqual(topology.num_clusters, 1)
        
        cell_radius = 1500
        num_clusters = 2
        topology = TopologySingleBaseStation(cell_radius, num_clusters)
        topology.calculate_coordinates()
        self.assertEqual(topology.num_clusters, 2)

        cell_radius = 1500
        num_clusters = 3
        with self.assertRaises(ValueError):
            topology = TopologySingleBaseStation(cell_radius, num_clusters)
            

    def test_intersite_distance(self):
        cell_radius = 1500
        num_clusters = 1
        topology = TopologySingleBaseStation(cell_radius, num_clusters)
        topology.calculate_coordinates()
        self.assertEqual(topology.intersite_distance, 3000)
        
        cell_radius = 1000
        num_clusters = 1
        topology = TopologySingleBaseStation(cell_radius, num_clusters)
        topology.calculate_coordinates()
        self.assertEqual(topology.intersite_distance, 2000)

        
    def test_cell_radius(self):
        cell_radius = 1500
        num_clusters = 1
        topology = TopologySingleBaseStation(cell_radius, num_clusters)
        topology.calculate_coordinates()
        self.assertEqual(topology.cell_radius, 1500)
        
        cell_radius = 1000
        num_clusters = 1
        topology = TopologySingleBaseStation(cell_radius, num_clusters)
        topology.calculate_coordinates()
        self.assertEqual(topology.cell_radius, 1000)
        
        
    def test_azimuth(self):
        cell_radius = 1500
        num_clusters = 1
        topology = TopologySingleBaseStation(cell_radius, num_clusters)
        topology.calculate_coordinates()
        self.assertEqual(topology.azimuth, 0)
        
        cell_radius = 1000
        num_clusters = 2
        topology = TopologySingleBaseStation(cell_radius, num_clusters)
        topology.calculate_coordinates()
        npt.assert_equal(topology.azimuth, np.array([0, 180]))
        

    def test_elevation(self):
        cell_radius = 1500
        num_clusters = 1
        topology = TopologySingleBaseStation(cell_radius, num_clusters)
        topology.calculate_coordinates()
        self.assertEqual(topology.elevation, -10)
        
        cell_radius = 1000
        num_clusters = 2
        topology = TopologySingleBaseStation(cell_radius, num_clusters)
        topology.calculate_coordinates()
        npt.assert_equal(topology.elevation, np.array([-10, -10]))
        
            
if __name__ == '__main__':
    unittest.main()        