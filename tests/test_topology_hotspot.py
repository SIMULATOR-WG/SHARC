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
        candidate_x = np.array([300])
        candidate_y = np.array([0])
        candidate_azimuth = np.array([-180])
        set_x = np.array([0, 200])
        set_y = np.array([0,   0])
        set_azimuth = np.array([0, -180])
        radius = 100
        self.assertFalse(self.topology.overlapping_hotspots(candidate_x, 
                                                            candidate_y, 
                                                            candidate_azimuth, 
                                                            set_x, 
                                                            set_y, 
                                                            set_azimuth,
                                                            radius))
        
        candidate_x = np.array([0])
        candidate_y = np.array([0])
        candidate_azimuth = np.array([0])
        set_x = np.array([0, 0])
        set_y = np.array([150, 400])
        set_azimuth = np.array([270, 270])
        radius = 100
        self.assertTrue(self.topology.overlapping_hotspots(candidate_x, 
                                                           candidate_y, 
                                                           candidate_azimuth, 
                                                           set_x, 
                                                           set_y, 
                                                           set_azimuth,
                                                           radius))
                
        candidate_x = np.array([0])
        candidate_y = np.array([0])
        candidate_azimuth = np.array([0])
        set_x = np.array([ -1, 101])
        set_y = np.array([  0,   0])
        set_azimuth = np.array([180, 0])
        radius = 100
        self.assertFalse(self.topology.overlapping_hotspots(candidate_x, 
                                                            candidate_y, 
                                                            candidate_azimuth, 
                                                            set_x, 
                                                            set_y, 
                                                            set_azimuth,
                                                            radius))
        
        candidate_x = np.array([1])
        candidate_y = np.array([0])
        candidate_azimuth = np.array([0])
        set_x = np.array([ 0])
        set_y = np.array([ 1])
        set_azimuth = np.array([90])
        radius = 100
        self.assertTrue(self.topology.overlapping_hotspots(candidate_x, 
                                                           candidate_y, 
                                                           candidate_azimuth, 
                                                           set_x, 
                                                           set_y, 
                                                           set_azimuth,
                                                           radius))
        
    
        
if __name__ == '__main__':
    unittest.main()        