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
        # Initialize variables for BS Cat-A mask (3.5 GHz)
        sta_type = StationType.IMT_BS
        p_tx = 46
        frequency = 3490
        bandwidth = 20
        spurious = -13
        self.mask_bs_a = SpectralMask3Gpp(sta_type, frequency, bandwidth, spurious)
        self.mask_bs_a.set_mask(power = p_tx)
        
        # Initialize variables for BS Cat-B mask (3.5 GHz)
        sta_type = StationType.IMT_BS
        p_tx = 46
        frequency = 3490
        bandwidth = 20
        spurious = -30
        self.mask_bs_b = SpectralMask3Gpp(sta_type, frequency, bandwidth, spurious)
        self.mask_bs_b.set_mask(power = p_tx)        
        
        # Initialize variables for UE mask (3.5 GHz)
        sta_type = StationType.IMT_UE
        p_tx = 22
        frequency = 3490
        bandwidth = 20
        spurious = -30
        self.mask_ue = SpectralMask3Gpp(sta_type, frequency, bandwidth, spurious)
        self.mask_ue.set_mask(power = p_tx)        
        
        
    def test_power_calc_bs_a(self):
        #######################################################################
        # BS Cat-A mask
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
        
    def test_power_calc_bs_b(self):        
        #######################################################################
        # BS Cat-B mask
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
        
    def test_power_calc_ue(self):        
        #######################################################################
        # UE mask
        #######################################################################
        fc = 3500.5
        bandwidth = 1
        poob = self.mask_ue.power_calc(fc, bandwidth)
        self.assertAlmostEqual(poob, -21 + 10*np.log10(1/0.03), delta = 1e-2)

        #######################################################################
        fc = 3503
        bandwidth = 4
        poob = self.mask_ue.power_calc(fc, bandwidth)
        self.assertAlmostEqual(poob, -10 + 10*np.log10(4), delta = 1e-2)
        
        ######################################################################
        fc = 3502.5
        bandwidth = 5
        poob = self.mask_ue.power_calc(fc, bandwidth)
        self.assertAlmostEqual(poob, 
                               10*np.log10(np.power(10, 0.1*(-21 + 10*np.log10(1/0.03))) + np.power(10, 0.1*(-10 + 10*np.log10(4)))), 
                               delta = 1e-2)
        
        #######################################################################
        fc = 3510
        bandwidth = 10
        poob = self.mask_ue.power_calc(fc, bandwidth)
        self.assertAlmostEqual(poob, 
                               -13 + 10*np.log10(10), 
                               delta = 1e-2)
        
        ######################################################################
        fc = 3520
        bandwidth = 10
        poob = self.mask_ue.power_calc(fc, bandwidth)
        self.assertAlmostEqual(poob, 
                               10*np.log10(np.power(10, 0.1*(-13 + 10*np.log10(5))) + np.power(10, 0.1*(-25 + 10*np.log10(5)))), 
                               delta = 1e-2)

        #######################################################################
        fc = 3525
        bandwidth = 10
        poob = self.mask_ue.power_calc(fc, bandwidth)
        self.assertAlmostEqual(poob, 
                               10*np.log10(np.power(10, 0.1*(-25 + 10*np.log10(5))) + np.power(10, 0.1*(-30 + 10*np.log10(5)))), 
                               delta = 1e-2)

        #######################################################################
        fc = 3600
        bandwidth = 50
        poob = self.mask_ue.power_calc(fc, bandwidth)
        self.assertAlmostEqual(poob, 
                               -30 + 10*np.log10(50),
                               delta = 1e-2)

        
if __name__ == '__main__':
    unittest.main()