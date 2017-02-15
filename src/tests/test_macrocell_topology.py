# -*- coding: utf-8 -*-
"""
Created on Tue Feb 14 13:05:42 2017

@author: edgar
"""

import unittest
import numpy as np

from macrocell_topology import MacrocellTopology

class MacrocellTopologyTest(unittest.TestCase):
    
    def setUp(self):
        intersite_distance = 1000
        num_clusters = 1
        self.topology = MacrocellTopology(intersite_distance, num_clusters)
        
    def test_intersite_distance(self):
        self.assertEqual(1000, self.topology.get_intersite_distance())
        self.topology.set_intersite_distance(600)
        self.assertEqual(600, self.topology.get_intersite_distance())
        self.topology.set_intersite_distance(1500)
        self.assertEqual(1500, self.topology.get_intersite_distance())

    def test_num_clusters(self):
        self.assertEqual(1, self.topology.get_num_clusters())
        # set to 7 clusters
        self.topology.set_num_clusters(7)
        self.assertEqual(7, self.topology.get_num_clusters())
        # set to 8 clusters, it raises an exception
        # NOTE: for some reason, running the assertRaises test crashes spyder
#        with self.assertRaises(ValueError):
#            self.topology.set_num_clusters(8)

    def test_coordinates(self):
        """
        TODO: test the case when number of clusters is 7 
        """
        (x_coord, y_coord) = self.topology.get_coordinates()
        # check the number of base stations
        self.assertEqual(19, len(x_coord))
        self.assertEqual(19, len(y_coord))
        # check coordinates
        x_ref = np.array([0, 1000, 500, -500, -1000, -500, 
                         500, 2000, 1500, 1000, 0, -1000, 
                         -1500, -2000, -1500, -1000, 0, 1000, 1500])
        y_ref = np.array([0, 0, 866.02, 866.02, 0, -866.02, 
                         -866.02, 0, 866.02, 1732.05, 1732.05, 1732.05, 
                         866.02, 0, -866.02, -1732.05, -1732.05, -1732.05, -866.02])
        self.assertTrue(np.all(np.isclose(x_coord, x_ref, atol=1e-2)))
        self.assertTrue(np.all(np.isclose(y_coord, y_ref, atol=1e-2)))
        
        # change intersite distance and check new coordinates
        self.topology.set_intersite_distance(500)
        (x_coord_h, y_coord_h) = self.topology.get_coordinates()
        # check the number of base stations
        self.assertEqual(19, len(x_coord_h))
        self.assertEqual(19, len(y_coord_h))
        # check coordinates
        x_ref_h = np.array([0, 500, 250, -250, -500, -250, 
                         250, 1000, 750, 500, 0, -500, 
                         -750, -1000, -750, -500, 0, 500, 750])
        y_ref_h = np.array([0, 0, 433.01, 433.01, 0, -433.01, 
                         -433.01, 0, 433.01, 866.02, 866.02, 866.02, 
                         433.01, 0, -433.01, -866.02, -866.02, -866.02, -433.01])
        self.assertTrue(np.all(np.isclose(x_coord_h, x_ref_h, atol=1e-2)))
        self.assertTrue(np.all(np.isclose(y_coord_h, y_ref_h, atol=1e-2)))       
            
if __name__ == '__main__':
    unittest.main()
        