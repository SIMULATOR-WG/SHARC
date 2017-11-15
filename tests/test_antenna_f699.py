#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 08:58:43 2017

@author: carlosrodriguez
"""

import unittest

from sharc.antenna.antenna_f699 import AntennaF699
from sharc.parameters.parameters_fs import ParametersFs

import numpy as np
import numpy.testing as npt

class AntennaF699Test(unittest.TestCase):
    
    def setUp(self):
        param = ParametersFs()
        param.antenna_gain = 40
        param.antenna_pattern = "ITU-R F.699-7"
        param.diameter = 0.3
        param.frequency = 27250

        self.antennad_lmdaless100 = AntennaF699(param)
        
        param.diameter = 2    
        self.antennad_lmdagreater100 = AntennaF699(param)
        
        
    def test_calculate_gain(self):
        psi = np.array([0,	0.1,	0.2,	0.5,	1,	2,	3,	4,	5,	10,	100])

        ref_gaind_lmdaless100 = np.array([40, 39.98,	39.92,	39.53,	38.14,	32.57,	23.53,	22.59,	20.17,	12.64,	-4.35])
        gaind_lmdaless100 = self.antennad_lmdaless100.calculate_gain(phi_vec=psi)
        npt.assert_allclose(gaind_lmdaless100, ref_gaind_lmdaless100, atol=1e-2)
    
        ref_gaind_lmdagreater100 = np.array([40,	39.17,	36.69,	35.88,	32,	24.47,	20.07,	16.94,	14.52,	7,	-10])
        gaind_lmdagreater100 = self.antennad_lmdagreater100.calculate_gain(phi_vec=psi) 
        npt.assert_allclose(gaind_lmdagreater100, ref_gaind_lmdagreater100, atol=1e-2)

        
if __name__ == '__main__':
    unittest.main()    