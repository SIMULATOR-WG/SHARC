# -*- coding: utf-8 -*-
"""
Created on Fri May  3 17:18:41 2019

@author: edgar
"""

import unittest
import numpy as np

from sharc.mask.spectral_mask_3gpp import SpectralMask3Gpp
from sharc.support.enumerations import StationType

class SpectalMask3GppTest(unittest.TestCase):
    
    def setUp(self):
        # Initialize variables for Cat-A mask (3.5 GHz)
        sta_type = StationType.IMT_BS
        p_tx = 46
        frequency = 3490
        bandwidth = 20
        spurious = -13
        self.mask_bs_a = SpectralMask3Gpp(sta_type, frequency, bandwidth, spurious)
        self.mask_bs_a.set_mask(power = p_tx)
        
        # Initialize variables for Cat-B mask (3.5 GHz)
        sta_type = StationType.IMT_BS
        p_tx = 46
        frequency = 3490
        bandwidth = 20
        spurious = -30
        self.mask_bs_b = SpectralMask3Gpp(sta_type, frequency, bandwidth, spurious)
        self.mask_bs_b.set_mask(power = p_tx)        
        
      
    def test_power_calc(self):
        #######################################################################
        # Cat-A mask
        #######################################################################
        fc = 3502.5
        bandwidth = 5
        poob = self.mask_bs_a.power_calc(fc, bandwidth)
        self.assertAlmostEqual(poob, 7.01, delta = 1e-2)

        #######################################################################
        fc = 3507.5
        bandwidth = 5
        poob = self.mask_bs_a.power_calc(fc, bandwidth)
        self.assertAlmostEqual(poob, 2.98, delta = 1e-2)
        
        #######################################################################
        fc = 3505
        bandwidth = 10
        poob = self.mask_bs_a.power_calc(fc, bandwidth)
        self.assertAlmostEqual(poob, 8.46, delta = 1e-2)
        
        #######################################################################
        fc = 3510
        bandwidth = 10
        poob = self.mask_bs_a.power_calc(fc, bandwidth)
        self.assertAlmostEqual(poob, 3.50, delta = 1e-2)
        
        #######################################################################
        fc = 3515
        bandwidth = 10
        poob = self.mask_bs_a.power_calc(fc, bandwidth)
        self.assertAlmostEqual(poob, -3, delta = 1e-2)

        #######################################################################
        fc = 3520
        bandwidth = 20
        poob = self.mask_bs_a.power_calc(fc, bandwidth)
        self.assertAlmostEqual(poob, 0.01, delta = 1e-2)
        
        #######################################################################
        # Cat-B mask
        #######################################################################
        fc = 3502.5
        bandwidth = 5
        poob = self.mask_bs_b.power_calc(fc, bandwidth)
        self.assertAlmostEqual(poob, 7.01, delta = 1e-2)

        #######################################################################
        fc = 3507.5
        bandwidth = 5
        poob = self.mask_bs_b.power_calc(fc, bandwidth)
        self.assertAlmostEqual(poob, 2.98, delta = 1e-2)
        
        #######################################################################
        fc = 3505
        bandwidth = 10
        poob = self.mask_bs_b.power_calc(fc, bandwidth)
        self.assertAlmostEqual(poob, 8.46, delta = 1e-2)
        
        #######################################################################
        fc = 3510
        bandwidth = 10
        poob = self.mask_bs_b.power_calc(fc, bandwidth)
        self.assertAlmostEqual(poob, 3, delta = 1e-2)
        
        #######################################################################
        fc = 3515
        bandwidth = 10
        poob = self.mask_bs_b.power_calc(fc, bandwidth)
        self.assertAlmostEqual(poob, -20, delta = 1e-2)

        #######################################################################
        fc = 3520
        bandwidth = 20
        poob = self.mask_bs_b.power_calc(fc, bandwidth)
        self.assertAlmostEqual(poob, -16.98, delta = 1e-2)        
        
        
if __name__ == '__main__':
    unittest.main()