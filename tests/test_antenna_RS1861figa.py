#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 10:31:53 2017

@author: carlosrodriguez
"""

import unittest

from sharc.antenna.antenna_rs1861_fig9a import AntennaRS1861FIG9a
from sharc.parameters.parameters_fss_ss import ParametersFssSs

import numpy as np
import numpy.testing as npt

class AntennaRS1861FIG9aTest(unittest.TestCase):
    
    def setUp(self):
        param = ParametersFssSs()
        param.antenna_pattern = "ITU-R RS.1861 figure a"
        param.antenna_gain = 50


        self.antenna = AntennaRS1861FIG9a(param)
        
        
        
    def test_calculate_gain(self):
        psi = np.array([0,	1, 2, 3, 4, 5, 6, 7, 8, 9, 	10,	100])

        ref_gain = np.array([50,42.01, 29.99,	15.09, 11.46,	7.82, 4.18, 0.54, -3.08, -6.72, -10, -10])
        gain = self.antenna.calculate_gain(phi_vec=psi)
        npt.assert_allclose(gain, ref_gain, atol=1e-2)
    
        
if __name__ == '__main__':
    unittest.main()   