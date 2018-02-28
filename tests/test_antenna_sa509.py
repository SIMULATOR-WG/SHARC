# -*- coding: utf-8 -*-
"""
Created on Thu Nov 16 17:47:45 2017

@author: Calil
"""

from sharc.antenna.antenna_sa509 import AntennaSA509
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

        self.antenna = AntennaSA509(self.par)

    def test_construction(self):
        self.assertEqual(self.antenna.diameter,10)
        self.assertEqual(self.antenna.efficiency,1)
        self.assertEqual(self.antenna.wavelength,1e-2)

        self.assertAlmostEqual(self.antenna.effective_area,78.539,delta=1e-2)

        self.assertAlmostEqual(self.antenna.g_0,69.943,delta=1e-2)
        self.assertAlmostEqual(self.antenna.phi_0,0.03464,delta=1e-4)

        self.assertAlmostEqual(self.antenna.phi_1,0.08944,delta=1e-4)
        self.assertAlmostEqual(self.antenna.phi_2,0.14531,delta=1e-4)

    def test_calculate_gain(self):
        phi = np.array([0.03464, 0.05, 0.1, 10, 25, 50, 100, 150])
        ref_gain = np.array([66.943, 63.693, 49.943, 4, -5.949, -13, -8, -13])

        gain = self.antenna.calculate_gain(off_axis_angle_vec=phi)
        npt.assert_allclose(gain, ref_gain, atol=1e-2)


if __name__ == '__main__':
    unittest.main()
