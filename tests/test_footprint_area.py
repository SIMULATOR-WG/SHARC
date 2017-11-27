# -*- coding: utf-8 -*-
"""
Created on Mon Nov 27 08:53:15 2017

@author: Calil
"""

import unittest
import numpy as np

from sharc.support.footprint_area import FootprintArea

class FootprintAreaTest(unittest.TestCase):
    
    def setUp(self):
        self.fa1 = FootprintArea(0,0,10)
        
    def test_construction(self):
        self.assertEqual(self.fa1.bore_lat_deg,0)
        self.assertEqual(self.fa1.bore_subsat_long_deg,0)
        self.assertEqual(self.fa1.beam_width,0)
        self.assertEqual(self.fa1.beta,0)
        self.assertEqual(self.fa1.bore_tilt,0)
        
if __name__ == '__main__':
    unittest.main()