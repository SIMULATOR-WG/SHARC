# -*- coding: utf-8 -*-
"""
Created on Sat Apr 15 15:36:22 2017

@author: Calil
"""

import unittest
import numpy as np

from sharc.antenna.antenna_beamforming_imt import AntennaBeamformingImt
from sharc.antenna.parameters.parameters_antenna_imt import ParametersAntennaImt

class AntennaBeamformingImtTest(unittest.TestCase):
    
    def setUp(self):
        #Array parameters
        self.param = ParametersAntennaImt()
        self.param.bs_element_max_g = 5
        self.param.bs_element_phi_3db = 80
        self.param.bs_element_theta_3db = 60
        self.param.bs_element_am = 30
        self.param.bs_element_sla_v = 30
        self.param.bs_n_rows = 16
        self.param.bs_n_columns = 16
        self.param.bs_element_horiz_spacing = 1
        self.param.bs_element_vert_spacing = 1
        self.param.ue_element_max_g = 10
        self.param.ue_element_phi_3db = 75
        self.param.ue_element_theta_3db = 65
        self.param.ue_element_am = 25
        self.param.ue_element_sla_v = 35
        self.param.ue_n_rows = 2
        self.param.ue_n_columns = 2
        self.param.ue_element_horiz_spacing = 0.5
        self.param.ue_element_vert_spacing = 0.5
        
        # Create antenna objects
        self.antenna1 = AntennaBeamformingImt(self.param,"BS")
        self.antenna2 = AntennaBeamformingImt(self.param,"UE")
        self.antenna2.gain = 15
        
    def test_gain(self):
        self.assertEqual(self.antenna1.gain,0)
        self.assertEqual(self.antenna2.gain,15)
        
    def test_station_type(self):
        self.assertTrue(self.antenna1.station_type == "BS")
        self.assertTrue(self.antenna2.station_type == "UE")
        
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
        
    def test_n_rows(self):
        self.assertEqual(self.antenna1.n_rows,16)
        self.assertEqual(self.antenna2.n_rows,2)
        
    def test_n_cols(self):
        self.assertEqual(self.antenna1.n_cols,16)
        self.assertEqual(self.antenna2.n_cols,2)
        
    def test_dh(self):
        self.assertEqual(self.antenna1.dh,1)
        self.assertEqual(self.antenna2.dh,0.5)
        
    def test_dv(self):
        self.assertEqual(self.antenna1.dv,1)
        self.assertEqual(self.antenna2.dv,0.5)
        
    def test_add(self):
        self.assertEqual(self.antenna1 + 2, 2)
    
    def test_radd(self):
        self.assertEqual(2 + self.antenna1, 2)
        self.assertEqual(self.antenna1 + self.antenna2, 15)
        
    def test_sub(self):
        self.assertEqual(self.antenna1 - 1, -1)
    
    def test_rsub(self):
        self.assertEqual(9 - self.antenna1, 9)
        self.assertEqual(self.antenna2 - self.antenna1, 15)
     
    def test_lt(self):
        self.assertTrue(self.antenna1 < 5)
        self.assertFalse(self.antenna1 < -3)
        self.assertTrue(self.antenna1 < self.antenna2)
        
    def test_le(self):
        self.assertTrue(self.antenna1 <= 5)
        self.assertTrue(self.antenna1 <= 0)
        self.assertFalse(self.antenna1 <= -2)
        self.assertTrue(self.antenna1 <= self.antenna2)

    def test_gt(self):
        self.assertFalse(self.antenna1 > 5)
        self.assertTrue(self.antenna1 > -2)
        self.assertFalse(self.antenna1 > self.antenna2)
        
    def test_ge(self):
        self.assertTrue(self.antenna1 >= 0)
        self.assertTrue(self.antenna1 >= -3)
        self.assertFalse(self.antenna1 >= 3)
        self.assertFalse(self.antenna1 >= self.antenna2)
        
    def test_eq(self):
        self.assertTrue(self.antenna1 == 0)
        self.assertFalse(self.antenna1 == 7)
        self.assertFalse(self.antenna1 == self.antenna2)
        
    def test_ne(self):
        self.assertTrue(self.antenna1 != 1)
        self.assertFalse(self.antenna1 != 0)
        self.assertTrue(self.antenna1 != self.antenna2)
        
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
        
    def test_element_pattern(self):
        # theta = 0 and phi = 90 result in maximum gain
        phi = 0
        theta = 90
        e_gain = self.antenna1.element_pattern(phi,theta)
        self.assertEqual(e_gain,5.0)
        self.assertEqual(e_gain,self.antenna1.g_max)
        
        # Check for gain update
        self.assertTrue(self.antenna1 == 5.0)
        
        phi = 80
        theta = 150
        e_gain = self.antenna1.element_pattern(phi,theta)
        self.assertEqual(e_gain,-19.0)
        self.assertTrue(self.antenna1 == -19.0)
        
        phi = 150
        theta = 210
        e_gain = self.antenna1.element_pattern(phi,theta)
        self.assertEqual(e_gain,-25.0)
        self.assertEqual(e_gain,self.antenna1.g_max - self.antenna1.am)
        self.assertTrue(self.antenna1 == -25.0)
        
    def test_super_position_vector(self):
        # Error margin
        eps = 1e-5
        
        # Test 1
        phi = 0
        theta = 0
        v_vec = self.antenna2.super_position_vector(phi, theta)
        expected_v_vec = np.array([[1.0, 1.0],[-1.0, -1.0]])
        self.assertTrue(np.allclose(np.real(v_vec),\
                                    np.real(expected_v_vec),rtol = eps))
        self.assertTrue(np.allclose(np.imag(v_vec),\
                                    np.imag(expected_v_vec),rtol = eps))
        
        # Test 2
        phi = 90
        theta = 90
        v_vec = self.antenna2.super_position_vector(phi, theta)
        expected_v_vec = np.array([[1.0, -1.0],[1.0, -1.0]])
        self.assertTrue(np.allclose(np.real(v_vec),\
                                    np.real(expected_v_vec),rtol = eps))
        self.assertTrue(np.allclose(np.imag(v_vec),\
                                    np.imag(expected_v_vec),rtol = eps))
        
        # Test 3
        phi = 45
        theta = 45
        v_vec = self.antenna2.super_position_vector(phi, theta)
        expected_v_vec = np.array([[1.0 + 0.0j, 0.0 + 1.0j],\
                    [-0.6056998+0.7956932j, -0.7956932-0.6056998j]])
        self.assertTrue(np.allclose(np.real(v_vec),\
                                    np.real(expected_v_vec),rtol = eps))
        self.assertTrue(np.allclose(np.imag(v_vec),\
                                    np.imag(expected_v_vec),rtol = eps))
        
        # Test 4
        phi = 60
        theta = 90
        v_vec = self.antenna2.super_position_vector(phi, theta)
        expected_v_vec = np.array([[1.0 + 0.0j, -0.912724 + 0.408576j],\
                    [1.0 + 0.0j, -0.912724 + 0.408576j]])
        self.assertTrue(np.allclose(np.real(v_vec),\
                                    np.real(expected_v_vec),rtol = eps))
        self.assertTrue(np.allclose(np.imag(v_vec),\
                                    np.imag(expected_v_vec),rtol = eps))
        
    def test_weight_vector(self):
        # Error margin
        eps = 1e-5
        
        # Test 1
        phi_scan = 0
        theta_tilt = 0
        w_vec = self.antenna2.weight_vector(phi_scan, theta_tilt)
        expected_w_vec = np.array([[0.5, 0.5],[0.5, 0.5]])
        self.assertTrue(np.allclose(np.real(w_vec),\
                                    np.real(expected_w_vec),rtol = eps))
        self.assertTrue(np.allclose(np.imag(w_vec),\
                                    np.imag(expected_w_vec),rtol = eps))
        
        # Test 2
        phi_scan = 90
        theta_tilt = 90
        w_vec = self.antenna2.weight_vector(phi_scan, theta_tilt)
        expected_w_vec = np.array([[0.5, 0.5],[-0.5, -0.5]])
        self.assertTrue(np.allclose(np.real(w_vec),\
                                    np.real(expected_w_vec),rtol = eps))
        self.assertTrue(np.allclose(np.imag(w_vec),\
                                    np.imag(expected_w_vec),rtol = eps))
        
        # Test 3
        phi_scan = 45
        theta_tilt = 45
        w_vec = self.antenna2.weight_vector(phi_scan, theta_tilt)
        expected_w_vec = np.array([[0.5 + 0.0j, 0.0 - 0.5j],\
                    [-0.3028499+0.3978466j, 0.3978466+0.3028499j]])
        self.assertTrue(np.allclose(np.real(w_vec),\
                                    np.real(expected_w_vec),rtol = eps))
        self.assertTrue(np.allclose(np.imag(w_vec),\
                                    np.imag(expected_w_vec),rtol = eps))
        
        # Test 4
        phi_scan = 0
        theta_tilt = 90
        w_vec = self.antenna2.weight_vector(phi_scan, theta_tilt)
        expected_w_vec = np.array([[0.5, 0.5],[-0.5, -0.5]])
        self.assertTrue(np.allclose(np.real(w_vec),\
                                    np.real(expected_w_vec),rtol = eps))
        self.assertTrue(np.allclose(np.imag(w_vec),\
                                    np.imag(expected_w_vec),rtol = eps))
        
        # Test 5
        phi_scan = 45
        theta_tilt = 30
        w_vec = self.antenna2.weight_vector(phi_scan, theta_tilt)
        expected_w_vec = np.array([[0.5 + 0.0j, -0.172870 - 0.469169j],\
                                   [0.0 + 0.5j,  0.469165 - 0.172870j]])
        self.assertTrue(np.allclose(np.real(w_vec),\
                                    np.real(expected_w_vec),rtol = eps))
        self.assertTrue(np.allclose(np.imag(w_vec),\
                                    np.imag(expected_w_vec),rtol = eps))
        
    def test_array_gain(self):
        # Error margin
        eps = 1e-5
        
        # Test 1
        v_vec = np.array([[0.5, 0.5],[0.5, 0.5]])
        w_vec = np.array([[0.5, 0.5],[0.5, 0.5]])
        ar_gain = self.antenna2.array_gain(v_vec,w_vec)
        self.assertEqual(ar_gain,0)
        
        # Test 2
        v_vec = np.array([[0.2  + 0.3j,  0.0 - 0.2j],\
                          [0.59 + 0.0j, -0.3 + 0.1j]])
        w_vec = np.array([[0.5  - 0.2j, -0.3 - 0.4j],\
                          [-0.1 + 1.0j,  0.0 + 0.5j]])
        ar_gain = self.antenna2.array_gain(v_vec,w_vec)
        self.assertAlmostEqual(ar_gain,-4.283598, delta = eps)
        
        # Test 3
        v_vec = np.array([[-0.7 + 0.2j, -5.0 + 0.3j],\
                          [ 0.2 + 1.0j,  0.0 - 0.9j]])
        w_vec = np.array([[ 0.6 + 0.6j,  2.0 + 0.1j],\
                          [-0.7 + 0.3j, -0.3 + 2.0j]])
        ar_gain = self.antenna2.array_gain(v_vec,w_vec)
        self.assertAlmostEqual(ar_gain,19.301796, delta = eps)
        
    def test_beam_gain(self):
        # Error margin
        eps = 1e-5
        
        # Test 1
        phi = 45
        theta = 45
        phi_scan = 45
        theta_tilt = 45
        beam_g = self.antenna2.beam_gain(phi,theta,phi_scan,theta_tilt)
        self.assertAlmostEqual(beam_g,1.594268,delta = eps)
        self.assertEqual(beam_g,self.antenna2.gain)
        
        # Test 1
        phi = 0
        theta = 60
        phi_scan = 0
        theta_tilt = 90
        beam_g = self.antenna2.beam_gain(phi,theta,phi_scan,theta_tilt)
        self.assertAlmostEqual(beam_g,10.454087,delta = eps)
        self.assertEqual(beam_g,self.antenna2.gain)
        
if __name__ == '__main__':
    unittest.main()