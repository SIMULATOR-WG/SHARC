# -*- coding: utf-8 -*-
"""
Created on Tue Feb 14 13:05:42 2017

@author: edgar
"""

import unittest
import numpy as np
import numpy.testing as npt

from sharc.topology_macrocell import TopologyMacrocell

class TopologyMacrocellTest(unittest.TestCase):
    
    def setUp(self):
        intersite_distance = 1000
        num_clusters = 1
        self.topology = TopologyMacrocell(intersite_distance, num_clusters)
        
    def test_intersite_distance(self):
        # check initial value
        self.assertEqual(self.topology.intersite_distance, 1000)
        self.assertAlmostEqual(self.topology.cell_radius, 666.667, places=2)
        
        # when intersite distance changes, cell radius also changes
        self.topology.intersite_distance = 600
        self.assertEqual(self.topology.intersite_distance, 600)
        self.assertEqual(self.topology.cell_radius, 400)
        
        self.topology.intersite_distance = 1500
        self.assertEqual(self.topology.intersite_distance, 1500)
        self.assertAlmostEqual(self.topology.cell_radius, 1000, places=2)

    def test_cell_radius(self):
        # check initial value
        self.assertEqual(self.topology.intersite_distance, 1000)
        self.assertAlmostEqual(self.topology.cell_radius, 666.667, places=2)
        
        # when cell radius changes, intersite distance also changes
        self.topology.cell_radius = 400
        self.assertEqual(self.topology.intersite_distance, 600)
        self.assertEqual(self.topology.cell_radius, 400)
        
        self.topology.cell_radius = 1000
        self.assertEqual(self.topology.intersite_distance, 1500)
        self.assertAlmostEqual(self.topology.cell_radius, 1000, places=2)
        
    def test_num_clusters(self):
        # check initial value
        self.assertEqual(self.topology.num_clusters, 1)
        # set to 7 clusters
        self.topology.num_clusters = 7
        self.assertEqual(self.topology.num_clusters, 7)
        # set to any other value raises an exception
        with self.assertRaises(ValueError):
            self.topology.num_clusters = 8
        # set to 5 clusters, it also raises an exception
        with self.assertRaises(ValueError):
            self.topology.num_clusters = 5            

    def test_coordinates(self):
        """
        TODO: test the case when number of clusters is 7 
        """
        # check the number of base stations
        self.assertEqual(len(self.topology.x), 19)
        self.assertEqual(len(self.topology.y), 19)
        # check coordinates
        x_ref = np.array([0, 1000, 500, -500, -1000, -500, 
                         500, 2000, 1500, 1000, 0, -1000, 
                         -1500, -2000, -1500, -1000, 0, 1000, 1500])
        y_ref = np.array([0, 0, 866.02, 866.02, 0, -866.02, 
                         -866.02, 0, 866.02, 1732.05, 1732.05, 1732.05, 
                         866.02, 0, -866.02, -1732.05, -1732.05, -1732.05, -866.02])
        npt.assert_allclose(self.topology.x, x_ref, atol=1e-2)
        npt.assert_allclose(self.topology.y, y_ref, atol=1e-2)
        
        # change intersite distance and check new cell radius and coordinates
        self.topology.intersite_distance = 500
        self.assertAlmostEqual(self.topology.cell_radius, 333.333, places=2)
        # check the number of base stations
        self.assertEqual(len(self.topology.x), 19)
        self.assertEqual(len(self.topology.y), 19)
        # check coordinates
        x_ref_h = np.array([0, 500, 250, -250, -500, -250, 
                         250, 1000, 750, 500, 0, -500, 
                         -750, -1000, -750, -500, 0, 500, 750])
        y_ref_h = np.array([0, 0, 433.01, 433.01, 0, -433.01, 
                         -433.01, 0, 433.01, 866.02, 866.02, 866.02, 
                         433.01, 0, -433.01, -866.02, -866.02, -866.02, -433.01])
        npt.assert_allclose(self.topology.x, x_ref_h, atol=1e-2)
        npt.assert_allclose(self.topology.y, y_ref_h, atol=1e-2)        
            
    def test_limits(self):
        self.assertAlmostEqual(self.topology.x_min, -2666.66, places=1)
        self.assertAlmostEqual(self.topology.x_max,  2666.66, places=1)
        self.assertAlmostEqual(self.topology.y_min, -2398.71, places=1)
        self.assertAlmostEqual(self.topology.y_max,  2398.71, places=1)
        
        # change inter site distance; scenario limits also have to change
        self.topology.intersite_distance = 500
        self.assertAlmostEqual(self.topology.cell_radius, 333.333, places=1)
        
        self.assertAlmostEqual(self.topology.x_min, -1333.33, places=1)
        self.assertAlmostEqual(self.topology.x_max,  1333.33, places=1)
        self.assertAlmostEqual(self.topology.y_min, -1199.35, places=1)
        self.assertAlmostEqual(self.topology.y_max,  1199.35, places=1)
        
if __name__ == '__main__':
    unittest.main()
        