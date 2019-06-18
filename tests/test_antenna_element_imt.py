# -*- coding: utf-8 -*-
"""
Created on Fri Apr 14 14:19:38 2017

@author: Calil
"""

import numpy as np
import unittest

from sharc.antenna.antenna_element_imt_m2101 import AntennaElementImtM2101
from sharc.parameters.parameters_antenna_imt import ParametersAntennaImt
from sharc.support.enumerations import StationType

class AntennaImtTest(unittest.TestCase):

    def setUp(self):
        #Element parameters
        self.param = ParametersAntennaImt()

        self.param.adjacent_antenna_model = "SINGLE_ELEMENT"
        self.param.bs_element_pattern = "M2101"
        self.param.ue_element_pattern = "M2101"
        self.param.bs_minimum_array_gain = -200
        self.param.ue_minimum_array_gain = -200
        self.param.normalization = False
        self.param.bs_downtilt = 0

        self.param.bs_normalization_file = None
        self.param.bs_element_max_g    = 5
        self.param.bs_element_phi_3db  = 80
        self.param.bs_element_theta_3db = 60
        self.param.bs_element_am       = 30
        self.param.bs_element_sla_v    = 30
        self.param.bs_n_rows           = 8
        self.param.bs_n_columns        = 8
        self.param.bs_element_horiz_spacing = 0.5
        self.param.bs_element_vert_spacing = 0.5
        self.param.bs_multiplication_factor = 12

        self.param.ue_normalization_file = None
        self.param.ue_element_max_g    = 10
        self.param.ue_element_phi_3db  = 75
        self.param.ue_element_theta_3db = 65
        self.param.ue_element_am       = 25
        self.param.ue_element_sla_v    = 35
        self.param.ue_n_rows           = 4
        self.param.ue_n_columns        = 4
        self.param.ue_element_horiz_spacing = 0.5
        self.param.ue_element_vert_spacing = 0.5
        self.param.ue_multiplication_factor = 12

        # Create antenna IMT objects
        par = self.param.get_antenna_parameters(StationType.IMT_BS)
        self.antenna1 = AntennaElementImtM2101(par)
        par = self.param.get_antenna_parameters(StationType.IMT_UE)
        self.antenna2 = AntennaElementImtM2101(par)

    def test_g_max(self):
        self.assertEqual(self.antenna1.g_max,5)
        self.assertEqual(self.antenna2.g_max,10)

    def test_phi_3db(self):
        self.assertEqual(self.antenna1.phi_3db,80)
        self.assertEqual(self.antenna2.phi_3db,75)

    def test_theta_3db(self):
        self.assertEqual(self.antenna1.theta_3db,60)
        self.assertEqual(self.antenna2.theta_3db,65)

    def test_am(self):
        self.assertEqual(self.antenna1.am,30)
        self.assertEqual(self.antenna2.am,25)

    def test_sla_v(self):
        self.assertEqual(self.antenna1.sla_v,30)
        self.assertEqual(self.antenna2.sla_v,35)

    def test_horizontal_pattern(self):
        # phi = 0 results in zero gain
        phi = 0
        h_att = self.antenna1.horizontal_pattern(phi)
        self.assertEqual(h_att,0.0)

        # phi = 120 implies horizontal gain of of -27 dB
        phi = 120
        h_att = self.antenna1.horizontal_pattern(phi)
        self.assertEqual(h_att,-27.0)

        # phi = 150, horizontal attenuation equals to the front-to-back ratio
        phi = 150
        h_att = self.antenna1.horizontal_pattern(phi)
        self.assertEqual(h_att,-30)
        self.assertEqual(h_att,-1.0*self.antenna1.am)

        # Test vector
        phi = np.array([0, 120, 150])
        h_att = self.antenna1.horizontal_pattern(phi)
        self.assertTrue(np.all(h_att == np.array([0.0,-27.0,-30.0])))

    def test_vertical_pattern(self):
        # theta = 90 results in zero gain
        theta = 90
        v_att = self.antenna1.vertical_pattern(theta)
        self.assertEqual(v_att,0.0)

        # theta = 180 implies vertical gain of -27 dB
        theta = 180
        v_att = self.antenna1.vertical_pattern(theta)
        self.assertEqual(v_att,-27.0)

        # theta = 210, vertical attenuation equals vertical sidelobe attenuation
        theta = 210
        v_att = self.antenna1.vertical_pattern(theta)
        self.assertEqual(v_att,-30)
        self.assertEqual(v_att,-1.0*self.antenna1.sla_v)

        # Test vector
        theta = np.array([90, 180, 210])
        v_att = self.antenna1.vertical_pattern(theta)
        self.assertTrue(np.all(v_att == np.array([0.0,-27.0,-30.0])))

    def test_element_pattern(self):
        # theta = 0 and phi = 90 result in maximum gain
        phi = 0
        theta = 90
        e_gain = self.antenna1.element_pattern(phi,theta)
        self.assertEqual(e_gain,5.0)
        self.assertEqual(e_gain,self.antenna1.g_max)

        phi = 80
        theta = 150
        e_gain = self.antenna1.element_pattern(phi,theta)
        self.assertEqual(e_gain,-19.0)

        phi = 150
        theta = 210
        e_gain = self.antenna1.element_pattern(phi,theta)
        self.assertEqual(e_gain,-25.0)
        self.assertEqual(e_gain,self.antenna1.g_max - self.antenna1.am)

        # Test vector
        phi = np.array([0,80,150])
        theta = np.array([90,150,210])
        e_gain = self.antenna1.element_pattern(phi,theta)
        self.assertTrue(np.all(e_gain == np.array([5.0,-19.0,-25.0])))

if __name__ == '__main__':
    unittest.main()
