#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 11:03:43 2017

@author: carlosrodriguez
"""

import unittest

from sharc.antenna.antenna_s465 import AntennaS465
from sharc.parameters.parameters_fss_es import ParametersFssEs

import numpy as np
import numpy.testing as npt

class AntennaS465Test(unittest.TestCase):
    
    def setUp(self):
        param = ParametersFssEs()
        param.antenna_gain = 40
        param.antenna_pattern = "ITU-R S465-6"
        param.diameter = 0.3
        param.frequency = 27250

        self.antennad_lmdaless50 = AntennaS465(param)
        
        param.diameter = 0.8   
        self.antennad_lmdagreater50 = AntennaS465(param)
        
        
    def test_calculate_gain(self):
        psi = np.array([0,	1,	2,	3,	4,	5,	6,	7,	8,	10,	40, 60])

        ref_gaind_lmdaless50 = np.array([40,	40,	40,	40,	16.94,	14.52,	12.54,	10.87,	9.42,	7,	-8.05,	-10])
        gaind_lmdaless50 = self.antennad_lmdaless50.calculate_gain(phi_vec=psi)
        npt.assert_allclose(gaind_lmdaless50, ref_gaind_lmdaless50, atol=1e-2)
    
        ref_gaind_lmdagreater50 = np.array([40,	40,	24.47,	20.07,	16.94,	14.52,	12.54,	10.87,	9.42,	7,	-8.05,	-10])
        gaind_lmdagreater50 = self.antennad_lmdagreater50.calculate_gain(phi_vec=psi) 
        npt.assert_allclose(gaind_lmdagreater50, ref_gaind_lmdagreater50, atol=1e-2)

        
if __name__ == '__main__':
    unittest.main()    