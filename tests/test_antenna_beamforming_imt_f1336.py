# -*- coding: utf-8 -*-
"""
Created on Fri May 31 15:57:20 2019

@author: edgar
"""

import unittest
import numpy as np
import numpy.testing as npt

from sharc.antenna.antenna_beamforming_imt import AntennaBeamformingImt
from sharc.parameters.parameters_antenna_imt import ParametersAntennaImt
from sharc.support.named_tuples import AntennaPar
from sharc.support.enumerations import StationType

class AntennaBeamformingImtF1336Test(unittest.TestCase):

    def setUp(self):
        #Array parameters
        self.param = ParametersAntennaImt()

        self.param.adjacent_antenna_model = "SINGLE_ELEMENT"
        self.param.bs_element_pattern = "F1336"
        self.param.bs_minimum_array_gain = -200
        self.param.normalization = False
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

        # Create antenna objects
        par = self.param.get_antenna_parameters(StationType.IMT_BS)
        self.antenna1 = AntennaBeamformingImt(par, 300, -10)


    def test_azimuth(self):
        self.assertEqual(self.antenna1.azimuth, 300)


    def test_elevation(self):
        self.assertEqual(self.antenna1.elevation, -10)


    def test_g_max(self):
        self.assertEqual(self.antenna1.element.g_max, 18)


    def test_phi_3db(self):
        self.assertEqual(self.antenna1.element.phi_3db, 65)


    def test_theta_3db(self):
        self.assertAlmostEqual(self.antenna1.element.theta_3db, 7.55, delta = 1e-2)


    def test_n_rows(self):
        self.assertEqual(self.antenna1.n_rows, 1)


    def test_n_cols(self):
        self.assertEqual(self.antenna1.n_cols, 1)


    def test_beams_list(self):
        self.assertEqual(len(self.antenna1.beams_list), 0)


    def test_w_vec_list(self):
        self.assertEqual(len(self.antenna1.w_vec_list), 0)


    def test_super_position_vector(self):
        phi = 0
        theta = 0
        v_vec = self.antenna1._super_position_vector(phi, theta)
        expected_v_vec = np.array([[1]])
        npt.assert_allclose(np.real(v_vec), 
                            np.real(expected_v_vec), 
                            atol = 1e-2)
        npt.assert_allclose(np.imag(v_vec), 
                            np.imag(expected_v_vec), 
                            atol = 1e-2)        

        phi = 90
        theta = 90
        v_vec = self.antenna1._super_position_vector(phi, theta)
        expected_v_vec = np.array([[1]])
        npt.assert_allclose(np.real(v_vec), 
                            np.real(expected_v_vec), 
                            atol = 1e-2)
        npt.assert_allclose(np.imag(v_vec), 
                            np.imag(expected_v_vec), 
                            atol = 1e-2)        
        phi = 45
        theta = 45
        v_vec = self.antenna1._super_position_vector(phi, theta)
        expected_v_vec = np.array([[1]])
        npt.assert_allclose(np.real(v_vec), 
                            np.real(expected_v_vec), 
                            atol = 1e-2)
        npt.assert_allclose(np.imag(v_vec), 
                            np.imag(expected_v_vec), 
                            atol = 1e-2)        

        phi = 60
        theta = 90
        v_vec = self.antenna1._super_position_vector(phi, theta)
        expected_v_vec = np.array([[1]])
        npt.assert_allclose(np.real(v_vec), 
                            np.real(expected_v_vec), 
                            atol = 1e-2)
        npt.assert_allclose(np.imag(v_vec), 
                            np.imag(expected_v_vec), 
                            atol = 1e-2)        


    def test_weight_vector(self):
        phi_scan = 0
        theta_tilt = 0
        w_vec = self.antenna1._weight_vector(phi_scan, theta_tilt)
        expected_w_vec = np.array([[1]])
        npt.assert_allclose(np.real(w_vec), 
                            np.real(expected_w_vec), 
                            atol = 1e-2)
        npt.assert_allclose(np.imag(w_vec), 
                            np.imag(expected_w_vec), 
                            atol = 1e-2)        

        phi_scan = 90
        theta_tilt = 90
        w_vec = self.antenna1._weight_vector(phi_scan, theta_tilt)
        expected_w_vec = np.array([[1]])
        npt.assert_allclose(np.real(w_vec), 
                            np.real(expected_w_vec), 
                            atol = 1e-2)
        npt.assert_allclose(np.imag(w_vec), 
                            np.imag(expected_w_vec), 
                            atol = 1e-2)        

        phi_scan = 45
        theta_tilt = 45
        w_vec = self.antenna1._weight_vector(phi_scan, theta_tilt)
        expected_w_vec = np.array([[1]])
        npt.assert_allclose(np.real(w_vec), 
                            np.real(expected_w_vec), 
                            atol = 1e-2)
        npt.assert_allclose(np.imag(w_vec), 
                            np.imag(expected_w_vec), 
                            atol = 1e-2)        

        phi_scan = 0
        theta_tilt = 90
        w_vec = self.antenna1._weight_vector(phi_scan, theta_tilt)
        expected_w_vec = np.array([[1]])
        npt.assert_allclose(np.real(w_vec), 
                            np.real(expected_w_vec), 
                            atol = 1e-2)
        npt.assert_allclose(np.imag(w_vec), 
                            np.imag(expected_w_vec), 
                            atol = 1e-2)        

        phi_scan = 45
        theta_tilt = 30
        w_vec = self.antenna1._weight_vector(phi_scan, theta_tilt)
        expected_w_vec = np.array([[1]])
        npt.assert_allclose(np.real(w_vec), 
                            np.real(expected_w_vec), 
                            atol = 1e-2)
        npt.assert_allclose(np.imag(w_vec), 
                            np.imag(expected_w_vec), 
                            atol = 1e-2)        


    def test_add_beam(self):
        eps = 1e-5
        par = self.param.get_antenna_parameters(StationType.IMT_BS)
        self.antenna1 = AntennaBeamformingImt(par, 0, 0)

        # Add first beam
        phi_scan = 45
        theta_tilt = 120
        self.antenna1.add_beam(phi_scan, theta_tilt)

        self.assertEqual(len(self.antenna1.beams_list), 1)
        self.assertEqual(len(self.antenna1.w_vec_list), 1)

        # Add second beam
        phi_scan = 90
        theta_tilt = 180.0
        self.antenna1.add_beam(phi_scan, theta_tilt)

        # Test beams_list
        self.assertEqual(len(self.antenna1.beams_list), 2)
        self.assertEqual(len(self.antenna1.w_vec_list), 2)

        # Test first beam
        self.assertAlmostEqual(self.antenna1.beams_list[0][0], 45, delta = eps)
        self.assertAlmostEqual(self.antenna1.beams_list[0][1], 30, delta = eps)

        w_vec = self.antenna1.w_vec_list[0]
        expected_w_vec = np.array([[1]])
        npt.assert_allclose(np.real(w_vec), 
                            np.real(expected_w_vec), 
                            atol = 1e-2)
        npt.assert_allclose(np.imag(w_vec), 
                            np.imag(expected_w_vec), 
                            atol = 1e-2) 

        # Test second beam
        self.assertEqual(self.antenna1.beams_list[1][0], 90)
        self.assertEqual(self.antenna1.beams_list[1][1], 90)

        w_vec = self.antenna1.w_vec_list[1]
        expected_w_vec = np.array([[1]])
        npt.assert_allclose(np.real(w_vec), 
                            np.real(expected_w_vec), 
                            atol = 1e-2)
        npt.assert_allclose(np.imag(w_vec), 
                            np.imag(expected_w_vec), 
                            atol = 1e-2) 

        # Reset beams and test
        self.antenna1.reset_beams()
        self.assertEqual(len(self.antenna1.beams_list), 0)
        self.assertEqual(len(self.antenna1.w_vec_list), 0)

    def test_beam_gain(self):
        # Error margin and antenna
        eps = 1e-2
        par = self.param.get_antenna_parameters(StationType.IMT_BS)
        self.antenna1 = AntennaBeamformingImt(par, 0, 0)

        # Test 1
        phi = 0
        theta = 90
        beam = 0
        phi_scan = 45
        theta_tilt = 135
        self.antenna1.add_beam(phi_scan, theta_tilt)
        beam_g = self.antenna1._beam_gain(phi, theta, beam)
        self.assertAlmostEqual(beam_g, 18, delta = eps)

        # Test 2
        phi = 30
        theta = 135
        beam = 1
        phi_scan = 45
        theta_tilt = 180
        self.antenna1.add_beam(phi_scan, theta_tilt)
        beam_g = self.antenna1._beam_gain(phi, theta, beam)
        self.assertAlmostEqual(beam_g, -1.48, delta = eps)

        # Test 3
        phi = 150
        theta = 180
        beam_g = self.antenna1._beam_gain(phi, theta)
        self.assertAlmostEqual(beam_g, -6.45, delta = eps)

    def test_calculate_gain(self):
        # Error margin and antenna
        eps = 1e-2
        par = self.param.get_antenna_parameters(StationType.IMT_BS)
        self.antenna1 = AntennaBeamformingImt(par,0,0)

        # Test 1
        phi_vec = np.array([0, 30])
        theta_vec = np.array([90, 100])
        gains = self.antenna1.calculate_gain(phi_vec = phi_vec, 
                                             theta_vec = theta_vec)
        npt.assert_allclose(gains,np.array([18, 4.52]), atol=eps)

        # Test 2
        phi = 0
        theta = 95
        phi_scan = 45
        theta_tilt = 180
        self.antenna1.add_beam(phi_scan, theta_tilt)
        beams_l = np.zeros_like(phi_vec, dtype=int)
        gains = self.antenna1.calculate_gain(phi_vec = phi, 
                                             theta_vec = theta,
                                             beams_l = beams_l)
        npt.assert_allclose(gains, np.array([12.74]), atol=eps)

        # Test 3
        phi = 5
        theta = 98
        gains = self.antenna1.calculate_gain(phi_vec = phi, 
                                             theta_vec = theta,
                                             co_channel = False)
        npt.assert_allclose(gains, np.array([6.81]), atol = eps)
        
        
if __name__ == '__main__':
    unittest.main()
