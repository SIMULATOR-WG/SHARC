#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 09:33:49 2017

@author: carlosrodriguez
"""

import unittest

from sharc.antenna.antenna_rs1813 import AntennaRS1813
from sharc.parameters.parameters_fss_ss import ParametersFssSs

import numpy as np
import numpy.testing as npt

class AntennaRS1813Test(unittest.TestCase):
    
    def setUp(self):
        param = ParametersFssSs()
        param.antenna_pattern = "ITU-R RS.1813-1"
        param.diameter = 0.6
        param.frequency = 23800

        self.antenna = AntennaRS1813(param)
        
        
        
    def test_calculate_gain(self):
        psi = np.array([0,	1, 2, 3, 4, 5, 6, 7, 8, 9, 	10,	100])

        ref_gain = np.array([61.27,	57.19,	44.96,	24.57,	9.56,	7.13,	5.15,	3.48,	2.03,	0.75,	-0.38,	-21.38])
        gain = self.antenna.calculate_gain(phi_vec=psi)
        npt.assert_allclose(gain, ref_gain, atol=1e-2)
    
        
if __name__ == '__main__':
    unittest.main()   