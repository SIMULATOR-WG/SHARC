# -*- coding: utf-8 -*-
"""
Created on Fri May 31 14:08:19 2019

@author: edgar
"""

import numpy as np
import unittest
import numpy.testing as npt

from sharc.antenna.antenna_element_imt_f1336 import AntennaElementImtF1336
from sharc.parameters.parameters_antenna_imt import ParametersAntennaImt
from sharc.support.enumerations import StationType

class AntennaElementImtF1336Test(unittest.TestCase):

    def setUp(self):
        #Element parameters
        self.param = ParametersAntennaImt()

        self.param.adjacent_antenna_model = "SINGLE_ELEMENT"
        self.param.bs_element_pattern = "F1336"
        self.param.bs_minimum_array_gain = -200
        self.param.bs_normalization = False
        self.param.bs_downtilt = 0

        self.param.bs_normalization_file = None
        self.param.bs_element_max_g    = 18
        self.param.bs_element_phi_3db  = 65
        self.param.bs_element_theta_3db = 0
        self.param.bs_n_rows           = 1
        self.param.bs_n_columns        = 1
        
        self.param.bs_element_am = 30
        self.param.bs_element_sla_v = 30
        self.param.bs_element_horiz_spacing = 0.5
        self.param.bs_element_vert_spacing = 0.5
        self.param.bs_multiplication_factor = 12

        # Create antenna IMT objects
        par = self.param.get_antenna_parameters(StationType.IMT_BS)
        self.antenna1 = AntennaElementImtF1336(par)

    def test_g_max(self):
        self.assertEqual(self.antenna1.g_max, 18)

    def test_phi_3db(self):
        self.assertEqual(self.antenna1.phi_3db, 65)

    def test_theta_3db(self):
        self.assertAlmostEqual(self.antenna1.theta_3db, 7.55, delta = 1e-2)

    def test_antenna_parameters(self):
        self.assertAlmostEqual(self.antenna1.k_a, 0.7, delta = 1e-2)
        self.assertAlmostEqual(self.antenna1.k_p, 0.7, delta = 1e-2)
        self.assertAlmostEqual(self.antenna1.k_h, 0.7, delta = 1e-2)
        self.assertAlmostEqual(self.antenna1.lambda_k_h, -1.87, delta = 1e-2)
        self.assertAlmostEqual(self.antenna1.k_v, 0.3, delta = 1e-2)
        self.assertAlmostEqual(self.antenna1.incline_factor, 18.45, delta = 1e-2)
        self.assertAlmostEqual(self.antenna1.x_k, 0.94, delta = 1e-2)
        self.assertAlmostEqual(self.antenna1.lambda_k_v, 4.60, delta = 1e-2)
        self.assertAlmostEqual(self.antenna1.g_hr_180, -24.45, delta = 1e-2)
        self.assertAlmostEqual(self.antenna1.g_hr_0, 0, delta = 1e-2)
                               
        
    def test_horizontal_pattern(self):
        phi = 0
        h_att = self.antenna1.horizontal_pattern(phi)
        self.assertAlmostEqual(h_att, 0, delta = 1e-2)

        phi = 30
        h_att = self.antenna1.horizontal_pattern(phi)
        self.assertAlmostEqual(h_att, -2.55, delta = 1e-2)

        phi = 150
        h_att = self.antenna1.horizontal_pattern(phi)
        self.assertAlmostEqual(h_att, -24.45, delta = 1e-2)

        # Test vector
        phi = np.array([0, 30, 150])
        h_att = self.antenna1.horizontal_pattern(phi)
        npt.assert_allclose(h_att, np.array([0, -2.55, -24.45]), atol = 1e-2)


    def test_vertical_pattern(self):
        theta = 90
        v_att = self.antenna1.vertical_pattern(theta)
        self.assertAlmostEqual(v_att, 0, delta = 1e-2)

        theta = 135
        v_att = self.antenna1.vertical_pattern(theta)
        self.assertAlmostEqual(v_att, -18.90, delta = 1e-2)

        theta = 180
        v_att = self.antenna1.vertical_pattern(theta)
        self.assertAlmostEqual(v_att, -24.45, delta = 1e-2)

        # Test vector
        theta = np.array([90, 135, 180])
        v_att = self.antenna1.vertical_pattern(theta)
        npt.assert_allclose(v_att, np.array([0, -18.90, -24.45]), atol = 1e-2)


    def test_element_pattern(self):
        phi = 0
        theta = 90
        e_gain = self.antenna1.element_pattern(phi, theta)
        self.assertEqual(e_gain, 18)
        self.assertEqual(e_gain, self.antenna1.g_max)

        phi = 30
        theta = 135
        e_gain = self.antenna1.element_pattern(phi, theta)
        self.assertAlmostEqual(e_gain, -1.48, delta = 1e-2)

        phi = 150
        theta = 180
        e_gain = self.antenna1.element_pattern(phi, theta)
        self.assertAlmostEqual(e_gain, -6.45, delta = 1e-2)

        # Test vector
        phi = np.array([0, 30, 150])
        theta = np.array([90, 135, 180])
        e_gain = self.antenna1.element_pattern(phi, theta)
        npt.assert_allclose(e_gain, np.array([18, -1.48, -6.45]), atol = 1e-2)

if __name__ == '__main__':
    unittest.main()
