# -*- coding: utf-8 -*-
"""
Created on Tue Feb 14 13:05:42 2017

@author: edgar
"""

import unittest
import numpy as np
import numpy.testing as npt

from sharc.topology.topology_macrocell import TopologyMacrocell

class TopologyMacrocellTest(unittest.TestCase):
    
    def setUp(self):
        pass
        
    def test_intersite_distance(self):
        intersite_distance = 1000
        num_clusters = 1
        topology = TopologyMacrocell(intersite_distance, num_clusters)
        self.assertEqual(topology.intersite_distance, 1000)
        self.assertAlmostEqual(topology.cell_radius, 666.667, places=2)
        
        # when intersite distance changes, cell radius also changes
        intersite_distance = 600
        num_clusters = 1
        topology = TopologyMacrocell(intersite_distance, num_clusters)
        self.assertEqual(topology.intersite_distance, 600)
        self.assertEqual(topology.cell_radius, 400)
        
        # let's change it one more time...
        intersite_distance = 1500
        num_clusters = 1
        topology = TopologyMacrocell(intersite_distance, num_clusters)
        self.assertEqual(topology.intersite_distance, 1500)
        self.assertAlmostEqual(topology.cell_radius, 1000, places=2)

    def test_num_clusters(self):
        # set to 1 cluster
        intersite_distance = 1000
        num_clusters = 1
        topology = TopologyMacrocell(intersite_distance, num_clusters)
        self.assertEqual(topology.num_clusters, 1)
        
        # set to 7 clusters
        intersite_distance = 1000
        num_clusters = 7
        topology = TopologyMacrocell(intersite_distance, num_clusters)
        self.assertEqual(topology.num_clusters, 7)
        
        # set to any other value raises an exception
        intersite_distance = 1000
        num_clusters = 3
        with self.assertRaises(ValueError):
            topology = TopologyMacrocell(intersite_distance, num_clusters)
        
        intersite_distance = 1000
        num_clusters = 8
        with self.assertRaises(ValueError):
            topology = TopologyMacrocell(intersite_distance, num_clusters)
            
    def test_coordinates(self):
        """
        TODO: test the case when number of clusters is 7 
        """
        intersite_distance = 1000
        num_clusters = 1
        topology = TopologyMacrocell(intersite_distance, num_clusters)
        
        num_sites = 19
        num_bs_per_site = 3
        num_bs = num_sites*num_bs_per_site
        
        # check the number of base stations
        self.assertEqual(len(topology.x), num_bs)
        self.assertEqual(len(topology.y), num_bs)
        self.assertEqual(len(topology.azimuth), num_bs)
        
        # check coordinates
        x_ref = np.repeat(np.array([0, 1000, 500, -500, -1000, -500, 
                         500, 2000, 1500, 1000, 0, -1000, 
                         -1500, -2000, -1500, -1000, 0, 1000, 1500]), num_bs_per_site)
        y_ref = np.repeat(np.array([0, 0, 866.02, 866.02, 0, -866.02, 
                         -866.02, 0, 866.02, 1732.05, 1732.05, 1732.05, 
                         866.02, 0, -866.02, -1732.05, -1732.05, -1732.05, -866.02]), num_bs_per_site)
        az_ref = np.tile([60, 180, 300], num_sites)
        
        npt.assert_allclose(topology.x, x_ref, atol=1e-2)
        npt.assert_allclose(topology.y, y_ref, atol=1e-2)
        npt.assert_allclose(topology.azimuth, az_ref, atol=1e-2)
        
        # change intersite distance and check new cell radius and coordinates
        intersite_distance = 500
        num_clusters = 1
        topology = TopologyMacrocell(intersite_distance, num_clusters)

        # check the number of base stations
        self.assertEqual(len(topology.x), num_bs)
        self.assertEqual(len(topology.y), num_bs)
        self.assertEqual(len(topology.azimuth), num_bs)

        # check coordinates
        x_ref = np.repeat(np.array([0, 500, 250, -250, -500, -250, 
                         250, 1000, 750, 500, 0, -500, 
                         -750, -1000, -750, -500, 0, 500, 750]), num_bs_per_site)
        y_ref = np.repeat(np.array([0, 0, 433.01, 433.01, 0, -433.01, 
                         -433.01, 0, 433.01, 866.02, 866.02, 866.02, 
                         433.01, 0, -433.01, -866.02, -866.02, -866.02, -433.01]), num_bs_per_site)
        az_ref = np.tile([60, 180, 300], num_sites)
        
        npt.assert_allclose(topology.x, x_ref, atol=1e-2)
        npt.assert_allclose(topology.y, y_ref, atol=1e-2)
        npt.assert_allclose(topology.azimuth, az_ref, atol=1e-2)      
            
    def test_limits(self):
        intersite_distance = 1000
        num_clusters = 1
        topology = TopologyMacrocell(intersite_distance, num_clusters)        
        
        self.assertAlmostEqual(topology.x_min, -2666.66, places=1)
        self.assertAlmostEqual(topology.x_max,  2666.66, places=1)
        self.assertAlmostEqual(topology.y_min, -2398.71, places=1)
        self.assertAlmostEqual(topology.y_max,  2398.71, places=1)
        
        # change inter site distance; scenario limits also have to change
        intersite_distance = 500
        num_clusters = 1
        topology = TopologyMacrocell(intersite_distance, num_clusters)  
        
        self.assertAlmostEqual(topology.x_min, -1333.33, places=1)
        self.assertAlmostEqual(topology.x_max,  1333.33, places=1)
        self.assertAlmostEqual(topology.y_min, -1199.35, places=1)
        self.assertAlmostEqual(topology.y_max,  1199.35, places=1)
        
        
if __name__ == '__main__':
    unittest.main()
        