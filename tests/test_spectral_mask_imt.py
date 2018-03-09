# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 11:56:10 2017

@author: Calil
"""

import unittest
import numpy as np

from sharc.spectral_mask_imt import SpectralMaskImt
from sharc.support.enumerations import StationType

class SpectalMaskImtTest(unittest.TestCase):
    
    def setUp(self):
        # Initialize variables
        sta_type = StationType.IMT_BS
        p_tx = 25.1
        freq = 43000
        band = 200
    
        # Create mask
        self.mask1 = SpectralMaskImt(sta_type,freq,band)
        self.mask1.set_mask(power = p_tx)
      
    def test_power_calc(self):
        # Test 1
        fc = 43000
        band = 200
        with self.assertWarns(RuntimeWarning):
            poob = self.mask1.power_calc(fc,band)
        self.assertAlmostEqual(poob,-np.inf,delta=1e-2)
        
        # Test 2
        fc = 43300
        band = 600
        poob = self.mask1.power_calc(fc,band)
        self.assertAlmostEqual(poob,11.8003,delta=1e-2)
        
        # Test 3
        fc = 43000
        band = 1200
        poob = self.mask1.power_calc(fc,band)
        self.assertAlmostEqual(poob,14.8106,delta=1e-2)
        
        # Test 4
        fc = 45000
        band = 1000
        poob = self.mask1.power_calc(fc, band)
        self.assertAlmostEqual(poob, 17, delta=1e-2)        
        
if __name__ == '__main__':
    unittest.main()