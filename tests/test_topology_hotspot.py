# -*- coding: utf-8 -*-
"""
Created on Tue May 23 16:43:09 2017

@author: edgar
"""

import unittest
import numpy as np
#import numpy.testing as npt

from sharc.parameters.parameters_hotspot import ParametersHotspot
from sharc.topology.topology_hotspot import TopologyHotspot

class TopologyHotspotTest(unittest.TestCase):
    
    
    def setUp(self):
        # For this test case, hotspot parameters are useless because we are
        # testing only the validation methods
        param = ParametersHotspot()
        param.num_hotspots_per_cell = 1
        param.max_dist_hotspot_ue = 100
        param.min_dist_hotspot_ue = 5
        param.min_dist_bs_hotspot = 105

        intersite_distance = 1000
        num_clusters = 1
        self.topology = TopologyHotspot(param, intersite_distance, num_clusters)
        
        
    def test_overlapping_hotspots(self):
        x = np.array([0, 200, 300])
        y = np.array([0,   0,   0])
        azimuth = np.array([0, -180, -180])
        radius = np.array([100, 100, 100])
        self.assertFalse(self.topology.overlapping_hotspots(x, y, azimuth, radius))
        
        x = np.array([0, 0, 0])
        y = np.array([0, 150, 400])
        azimuth = np.array([90, 270, 270])
        radius = np.array([100, 100, 100])
        self.assertTrue(self.topology.overlapping_hotspots(x, y, azimuth, radius))        
                
        x = np.array([ 0, -1, 101])
        y = np.array([ 0,  0,   0])
        azimuth = np.array([0, 180, 0])
        radius = np.array([100, 100, 100])
        self.assertFalse(self.topology.overlapping_hotspots(x, y, azimuth, radius))       
        
        x = np.array([ 1, 0])
        y = np.array([ 0, 1])
        azimuth = np.array([0, 90])
        radius = np.array([100, 100])
        self.assertTrue(self.topology.overlapping_hotspots(x, y, azimuth, radius))         
        
    
        
if __name__ == '__main__':
    unittest.main()        