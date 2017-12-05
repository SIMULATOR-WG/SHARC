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
        df = -10
        band = 510
        power = self.mask1.power_calc(df,band)
        self.assertAlmostEqual(power,35.1322,delta=1e+2)
        
        # Test 2
        df = 410
        band = 10
        power = self.mask1.power_calc(df,band)
        self.assertAlmostEqual(power,-2.9999,delta=1e+2)
        
        # Test 3
        df = -20
        band = 10
        power = self.mask1.power_calc(df,band)
        self.assertAlmostEqual(power,35.1000,delta=1e+2)
        
if __name__ == '__main__':
    unittest.main()