#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 11:12:59 2017

@author: carlosrodriguez
"""

import unittest

from sharc.antenna.antenna_s509 import AntennaS509
from sharc.parameters.parameters_fss_es import ParametersFssEs

import numpy as np
import numpy.testing as npt

class AntennaS465Test(unittest.TestCase):
    
    def setUp(self):
        param = ParametersFssEs()
        param.antenna_gain = 40
        param.antenna_pattern = "ITU-R SA.509-3"
        param.diameter = 2
        param.frequency = 23800

        self.antenna = AntennaS509(param)
        
        
        
    def test_calculate_gain(self):
        psi = np.array([0,	1,	2,	3,	4,	5,	6,	7,	8,	9,	10, 100])

        ref_gain = np.array([40,	20,	20,	17.07,	13.94,	11.52,	9.54,	7.87,	6.42,	5.14,	4,	-8])
        gain = self.antenna.calculate_gain(phi_vec=psi)
        npt.assert_allclose(gain, ref_gain, atol=1e-2)
    
 
        
if __name__ == '__main__':
    unittest.main()   