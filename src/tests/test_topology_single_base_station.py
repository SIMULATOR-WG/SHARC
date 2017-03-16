# -*- coding: utf-8 -*-
"""
Created on Mon Mar 13 17:21:44 2017

@author: edgar
"""

import unittest

from topology_single_base_station import TopologySingleBaseStation

class TopologySingleBaseStationTest(unittest.TestCase):
    
    def setUp(self):
        cell_radius = 1500
        self.topology = TopologySingleBaseStation(cell_radius)
        
    def test_cell_radius(self):
        # check initial value
        self.assertEqual(self.topology.cell_radius, 1500)
        # change cell radius and check it
        self.topology.cell_radius = 700
        self.assertEqual(self.topology.cell_radius, 700)

        self.topology.cell_radius = 400
        self.assertEqual(self.topology.cell_radius, 400)
        
    def test_coordinates(self):
        self.assertEqual(self.topology.x, 0)
        self.assertEqual(self.topology.y, 0)
        
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
        
    def test_num_clusters(self):
        with self.assertRaises(NotImplementedError):
            self.topology.num_clusters = 1
        with self.assertRaises(NotImplementedError):
            self.topology.num_clusters = 2            
        with self.assertRaises(NotImplementedError):
            self.topology.num_clusters = 7
        
    def test_intersite_distance(self):
        with self.assertRaises(NotImplementedError):
            self.topology.intersite_distance = 100
        with self.assertRaises(NotImplementedError):
            self.topology.intersite_distance = 500        
        with self.assertRaises(NotImplementedError):
            self.topology.intersite_distance = 1000  
            
if __name__ == '__main__':
    unittest.main()        