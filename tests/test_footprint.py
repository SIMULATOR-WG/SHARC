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
        self.fa1 = Footprint(0,0,0.1)
        self.fa2 = Footprint(0,0,0.325)
        
    def test_construction(self):
        self.assertEqual(self.fa1.bore_lat_deg,0)
        self.assertEqual(self.fa1.bore_subsat_long_deg,0)
        self.assertEqual(self.fa1.beam_width_deg,0.1)
        self.assertEqual(self.fa1.bore_lat_rad,0)
        self.assertEqual(self.fa1.bore_subsat_long_rad,0)
        self.assertEqual(self.fa1.beam_width_rad,np.pi/1800)
        self.assertEqual(self.fa1.beta,0)
        self.assertEqual(self.fa1.bore_tilt,0)
        
    def test_calc_footprint(self):
        fp_long, fp_lat = self.fa1.calc_footprint(4)
        npt.assert_allclose(fp_long,np.array([-0.398, 0.543,  0.146, 0.398]),atol=1e-2)
        npt.assert_allclose(fp_lat,np.array([-0.398, -0.146, -0.543,  0.398]),atol=1e-2)
        
    def test_calc_area(self):
        a = self.fa2.calc_area(500)
        self.assertAlmostEqual(a,129464,delta=150)
        
if __name__ == '__main__':
    unittest.main()