# -*- coding: utf-8 -*-
"""
Created on Mon Mar 13 17:21:44 2017

@author: edgar
"""

import unittest
import numpy.testing as npt

from topology_single_base_station import TopologySingleBaseStation

class TopologySingleBaseStationTest(unittest.TestCase):
    
    def setUp(self):
        cell_radius = 1500
        num_clusters = 1
        self.topology = TopologySingleBaseStation(cell_radius, num_clusters)
        
    def test_coordinates(self):
        self.assertEqual(self.topology.x, [0])
        self.assertEqual(self.topology.y, [0])
        # with a single cluster there is no change in coordinates
        self.topology.intersite_distance = 1000
        self.assertEqual(self.topology.x, [0])
        self.assertEqual(self.topology.y, [0])
        self.topology.cell_radius = 1000
        self.assertEqual(self.topology.x, [0])
        self.assertEqual(self.topology.y, [0])
        
        # now we change to 2 clusters
        self.topology.cell_radius = 1500
        self.topology.num_clusters = 2
        npt.assert_array_equal(self.topology.x, [-1500, 1500])
        npt.assert_array_equal(self.topology.y, [0, 0])
        # 
        self.topology.intersite_distance = 1000
        npt.assert_array_equal(self.topology.x, [-500, 500])
        npt.assert_array_equal(self.topology.y, [0, 0])
        self.topology.cell_radius = 500
        npt.assert_array_equal(self.topology.x, [-500, 500])
        npt.assert_array_equal(self.topology.y, [0, 0])
        
    def test_limits(self):
        # check initial values
        self.assertEqual(self.topology.x_min, -1500)
        self.assertEqual(self.topology.x_max,  1500)
        self.assertEqual(self.topology.y_min, -1500)
        self.assertEqual(self.topology.y_max,  1500)
        # change cell radius and check it
        self.topology.cell_radius = 1000
        self.assertEqual(self.topology.x_min, -1000)
        self.assertEqual(self.topology.x_max,  1000)
        self.assertEqual(self.topology.y_min, -1000)
        self.assertEqual(self.topology.y_max,  1000)
        
        # change to 2 clusters
        self.topology.cell_radius = 1500
        self.topology.num_clusters = 2        
        self.assertEqual(self.topology.x_min, -3000)
        self.assertEqual(self.topology.x_max,  3000)
        self.assertEqual(self.topology.y_min, -1500)
        self.assertEqual(self.topology.y_max,  1500)
        # change cell radius and check it
        self.topology.cell_radius = 1000
        self.assertEqual(self.topology.x_min, -2000)
        self.assertEqual(self.topology.x_max,  2000)
        self.assertEqual(self.topology.y_min, -1000)
        self.assertEqual(self.topology.y_max,  1000)        
        
    def test_num_clusters(self):
        self.assertEqual(self.topology.num_clusters, 1)
        self.topology.num_clusters = 2
        self.assertEqual(self.topology.num_clusters, 2)
        with self.assertRaises(ValueError):
            self.topology.num_clusters = 3
        
    def test_intersite_distance(self):
        self.assertEqual(self.topology.intersite_distance, 3000)
        
        self.topology.intersite_distance = 1000
        self.assertEqual(self.topology.intersite_distance, 1000)
        self.assertEqual(self.topology.cell_radius, 500)
        
        self.topology.intersite_distance = 2000
        self.assertEqual(self.topology.intersite_distance, 2000)
        self.assertEqual(self.topology.cell_radius, 1000)      
        
    def test_cell_radius(self):
        # check initial value
        self.assertEqual(self.topology.cell_radius, 1500)
        # change cell radius and check it
        self.topology.cell_radius = 1000
        self.assertEqual(self.topology.intersite_distance, 2000)
        self.assertEqual(self.topology.cell_radius, 1000)
        # change cell radius and check it
        self.topology.cell_radius = 2000
        self.assertEqual(self.topology.intersite_distance, 4000)
        self.assertEqual(self.topology.cell_radius, 2000)         
            
if __name__ == '__main__':
    unittest.main()        