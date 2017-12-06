# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 11:56:10 2017

@author: Calil
"""

import unittest

from sharc.spectral_mask_imt import SpectralMaskImt
from sharc.support.enumerations import StationType

class SpectalMaskImtTest(unittest.TestCase):
    
    def setUp(self):
        # Initialize variables
        sta_type = StationType.IMT_BS
        p_tx = 25.1
        freq = 43000
    
        # Create mask
        self.mask1 = SpectralMaskImt(sta_type,p_tx,freq)
      
    def test_power_calc(self):
        # Test 1
        fc = 43000
        band = 200
        pib, poob = self.mask1.power_calc(fc,band)
        self.assertAlmostEqual(pib,25.1,delta=1e-2)
        self.assertAlmostEqual(poob,-500,delta=1e-2)
        
        # Test 2
        fc = 43300
        band = 600
        pib, poob = self.mask1.power_calc(fc,band)
        self.assertAlmostEqual(pib,22.1,delta=5e-2)
        self.assertAlmostEqual(poob,11.8003,delta=1e-2)
        
        # Test 3
        fc = 43000
        band = 1200
        pib, poob = self.mask1.power_calc(fc,band)
        self.assertAlmostEqual(pib,25.1,delta=5e-2)
        self.assertAlmostEqual(poob,14.8106,delta=1e-2)
        
if __name__ == '__main__':
    unittest.main()