# -*- coding: utf-8 -*-
"""
Created on Mon Nov 27 08:53:15 2017

@author: Calil
"""

import unittest
import numpy as np
import numpy.testing as npt

from sharc.support.footprint import Footprint

class FootprintAreaTest(unittest.TestCase):
    
    def setUp(self):
        self.fa1 = Footprint(0.1,bore_lat_deg=0,bore_subsat_long_deg=0.0)
        self.fa2 = Footprint(0.325,bore_lat_deg = 0)
        self.fa3 = Footprint(0.325,elevation_deg = 20)
        
    def test_construction(self):
        self.assertEqual(self.fa1.bore_lat_deg,0)
        self.assertEqual(self.fa1.bore_subsat_long_deg,0)
        self.assertEqual(self.fa1.beam_width_deg,0.1)
        self.assertEqual(self.fa1.bore_lat_rad,0)
        self.assertEqual(self.fa1.bore_subsat_long_rad,0)
        self.assertEqual(self.fa1.beam_width_rad,np.pi/1800)
        self.assertEqual(self.fa1.beta,0)
        self.assertEqual(self.fa1.bore_tilt,0)
        
        self.assertEqual(self.fa2.bore_lat_deg,0)
        self.assertEqual(self.fa2.bore_subsat_long_deg,0)
        self.assertEqual(self.fa2.bore_lat_rad,0)
        self.assertEqual(self.fa2.bore_subsat_long_rad,0)
        
        self.assertEqual(self.fa3.bore_lat_deg,0)
        self.assertAlmostEqual(self.fa3.bore_subsat_long_deg,61.84,delta=0.01)
        
    def test_set_elevation(self):
        self.fa2.set_elevation(20)
        self.assertEqual(self.fa2.bore_lat_deg,0)
        self.assertAlmostEqual(self.fa2.bore_subsat_long_deg,61.84,delta=0.01)
        
    def test_calc_footprint(self):
        fp_long, fp_lat = self.fa1.calc_footprint(4)
        npt.assert_allclose(fp_long,np.array([0.0, 0.487,  -0.487, 0.0]),atol=1e-2)
        npt.assert_allclose(fp_lat,np.array([-0.562,  0.281,  0.281,  -0.562]),atol=1e-2)
        
    def test_calc_area(self):
        a1 = self.fa2.calc_area(1000)
        self.assertAlmostEqual(a1,130000,delta=200)
        a2 = self.fa3.calc_area(1000)
        self.assertAlmostEqual(a2,486300,delta=200)
        
if __name__ == '__main__':
    unittest.main()