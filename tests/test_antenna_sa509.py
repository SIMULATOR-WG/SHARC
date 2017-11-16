# -*- coding: utf-8 -*-
"""
Created on Thu Nov 16 17:47:45 2017

@author: Calil
"""

from sharc.antenna.antenna_sa509 import AntenaSA509
from sharc.parameters.parameters_ras import ParametersRas

import numpy as np
import unittest
import numpy.testing as npt

class AntennaSA509Test(unittest.TestCase):
    
    def setUp(self):
        self.par = ParametersRas();
        self.par.diameter = 10
        self.par.antenna_efficiency = 1
        self.par.frequency = 30000
        self.par.SPEED_OF_LIGHT = 3e8
        
        self.antenna = AntenaSA509(self.par)
        
    def test_construction(self):
        self.assertEqual(self.antenna.diameter,10)
        self.assertEqual(self.antenna.efficiency,1)
        self.assertEqual(self.antenna.wavelength,1e-2)
        
        self.assertAlmostEqual(self.antenna.g_0,69.943,delta=1e-2)
        self.assertAlmostEqual(self.antenna.phi_0,0.03464,delta=1e-4)
        
        self.assertAlmostEqual(self.antenna.phi_1,0.08944,delta=1e-4)
        self.assertAlmostEqual(self.antenna.phi_2,0.14531,delta=1e-4)
        
        
if __name__ == '__main__':
    unittest.main()