# -*- coding: utf-8 -*-
"""
Created on Tue May 23 16:43:09 2017

@author: edgar
"""

import unittest
import numpy as np
#import numpy.testing as npt

from sharc.parameters.parameters_hotspot import ParametersHotspot
from sharc.topology.topology_macro_hotspot import TopologyMacroHotspot

class TopologyMacroHotspotTest(unittest.TestCase):
    
    def setUp(self):
        param = ParametersHotspot()
        param.num_hotspots_per_cell = 2
        param.num_stations_per_hotspots = 4
        param.ue_hotspot_dropping_ratio = 0.67
        param.ue_outdoor_ratio = 0.8
        param.max_dist_station_hotspot = 50
        param.max_dist_ue_hotspot = 70
        param.min_dist_station_station = 20
        param.min_dist_station_ue = 5
        param.min_dist_bs_hotspot = 105
        param.min_dist_bs_ue = 35
        param.min_dist_hotspots = 2*param.max_dist_station_hotspot

        intersite_distance = 1000
        num_clusters = 1
        self.topology = TopologyMacroHotspot(param, intersite_distance, num_clusters)
        
    def test_validade_min_dist_hotspots(self):
        hotspot_center_x = np.array([0, 100, 300])
        hotspot_center_y = np.array([0, 0, 0])
        min_dist_hotspots = 100
        self.assertTrue(self.topology.validade_min_dist_hotspots(hotspot_center_x, 
                                                                 hotspot_center_y, 
                                                                 min_dist_hotspots))
        
        hotspot_center_x = np.array([-100, 100, 400])
        hotspot_center_y = np.array([0, 0, 0])
        min_dist_hotspots = 150
        self.assertTrue(self.topology.validade_min_dist_hotspots(hotspot_center_x, 
                                                                 hotspot_center_y, 
                                                                 min_dist_hotspots))   
        
        hotspot_center_x = np.array([-100, 100, 400])
        hotspot_center_y = np.array([-50, 25, 200])
        min_dist_hotspots = 150
        self.assertTrue(self.topology.validade_min_dist_hotspots(hotspot_center_x, 
                                                                 hotspot_center_y, 
                                                                 min_dist_hotspots))         
        
        hotspot_center_x = np.array([-100, 100, 400])
        hotspot_center_y = np.array([0, 0, 0])
        min_dist_hotspots = 250
        self.assertFalse(self.topology.validade_min_dist_hotspots(hotspot_center_x, 
                                                                 hotspot_center_y, 
                                                                 min_dist_hotspots))         
        
        
    def test_validade_min_dist_bs_hotspot(self):
        hotspot_center_x = np.array([0, 100, 300])
        hotspot_center_y = np.array([0, 0, 0])        
        macrocell_x = np.array([-300, -100])
        macrocell_y = np.array([0, 0])   
        min_dist_bs_hotspot = 100
        self.assertTrue(self.topology.validade_min_dist_bs_hotspot(hotspot_center_x, 
                                                                   hotspot_center_y, 
                                                                   macrocell_x, 
                                                                   macrocell_y, 
                                                                   min_dist_bs_hotspot))
        
        hotspot_center_x = np.array([0, 100, 300])
        hotspot_center_y = np.array([0, 100, 300])        
        macrocell_x = np.array([-300, -100, 50])
        macrocell_y = np.array([0, 0, 50])   
        min_dist_bs_hotspot = 100
        self.assertFalse(self.topology.validade_min_dist_bs_hotspot(hotspot_center_x, 
                                                                    hotspot_center_y, 
                                                                    macrocell_x, 
                                                                    macrocell_y, 
                                                                    min_dist_bs_hotspot))
        
        hotspot_center_x = np.array([0, 100, 300, 1000])
        hotspot_center_y = np.array([0, 100, 300, 500])        
        macrocell_x = np.array([-300, -100, 50, 400])
        macrocell_y = np.array([0, 0, 50, 400])   
        min_dist_bs_hotspot = 50
        self.assertTrue(self.topology.validade_min_dist_bs_hotspot(hotspot_center_x, 
                                                                   hotspot_center_y, 
                                                                   macrocell_x, 
                                                                   macrocell_y, 
                                                                   min_dist_bs_hotspot))          
        
    def test_validade_min_dist_stations(self):
        station_x = np.array([0, 100, 300])
        station_y = np.array([0, 0, 0])
        min_dist_stations = 20
        self.assertTrue(self.topology.validade_min_dist_hotspots(station_x, 
                                                                 station_y, 
                                                                 min_dist_stations))
        
        station_x = np.array([0, 0, 0])
        station_y = np.array([0, 50, 80])
        min_dist_stations = 30
        self.assertTrue(self.topology.validade_min_dist_hotspots(station_x, 
                                                                 station_y, 
                                                                 min_dist_stations))
        
        station_x = np.array([0, 80, 300])
        station_y = np.array([0, 0, 0])
        min_dist_stations = 100
        self.assertFalse(self.topology.validade_min_dist_hotspots(station_x, 
                                                                  station_y, 
                                                                  min_dist_stations))
        
        station_x = np.array([0, 20, 100])
        station_y = np.array([0, 20, 0])
        min_dist_stations = 50
        self.assertFalse(self.topology.validade_min_dist_hotspots(station_x, 
                                                                  station_y, 
                                                                  min_dist_stations))
        
        
    def test_calculate_coordinates(self):
        self.topology.calculate_coordinates()
        
if __name__ == '__main__':
    unittest.main()        