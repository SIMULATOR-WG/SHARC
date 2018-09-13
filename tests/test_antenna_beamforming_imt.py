# -*- coding: utf-8 -*-
"""
Created on Sat Apr 15 15:36:22 2017

@author: Calil
"""

import unittest
import numpy as np
import numpy.testing as npt

from sharc.antenna.antenna_beamforming_imt import AntennaBeamformingImt
from sharc.parameters.parameters_antenna_imt import ParametersAntennaImt
from sharc.support.named_tuples import AntennaPar

class AntennaBeamformingImtTest(unittest.TestCase):

    def setUp(self):
        #Array parameters
        self.param = ParametersAntennaImt()

        self.param.normalization = False
        self.param.bs_normalization_file = None
        self.param.bs_element_pattern = "M2101"
        self.param.bs_downtilt_deg = 0

        self.param.bs_tx_element_max_g    = 5
        self.param.bs_tx_element_phi_deg_3db  = 80
        self.param.bs_tx_element_theta_deg_3db = 60
        self.param.bs_tx_element_am       = 30
        self.param.bs_tx_element_sla_v    = 30
        self.param.bs_tx_n_rows           = 16
        self.param.bs_tx_n_columns        = 16
        self.param.bs_tx_element_horiz_spacing = 1
        self.param.bs_tx_element_vert_spacing = 1

        self.param.bs_rx_element_max_g    = 5
        self.param.bs_rx_element_phi_deg_3db  = 80
        self.param.bs_rx_element_theta_deg_3db = 60
        self.param.bs_rx_element_am       = 30
        self.param.bs_rx_element_sla_v    = 30
        self.param.bs_rx_n_rows           = 16
        self.param.bs_rx_n_columns        = 16
        self.param.bs_rx_element_horiz_spacing = 1
        self.param.bs_rx_element_vert_spacing = 1

        self.param.ue_element_pattern = "M2101"
        self.param.ue_normalization_file = None

        self.param.ue_tx_element_max_g    = 10
        self.param.ue_tx_element_phi_deg_3db  = 75
        self.param.ue_tx_element_theta_deg_3db = 65
        self.param.ue_tx_element_am       = 25
        self.param.ue_tx_element_sla_v    = 35
        self.param.ue_tx_n_rows           = 2
        self.param.ue_tx_n_columns        = 2
        self.param.ue_tx_element_horiz_spacing = 0.5
        self.param.ue_tx_element_vert_spacing = 0.5

        self.param.ue_rx_element_max_g    = 10
        self.param.ue_rx_element_phi_deg_3db  = 75
        self.param.ue_rx_element_theta_deg_3db = 65
        self.param.ue_rx_element_am       = 25
        self.param.ue_rx_element_sla_v    = 35
        self.param.ue_rx_n_rows           = 2
        self.param.ue_rx_n_columns        = 2
        self.param.ue_rx_element_horiz_spacing = 0.5
        self.param.ue_rx_element_vert_spacing = 0.5

        # Create antenna objects
        par = self.param.get_antenna_parameters("BS","RX")
        self.antenna1 = AntennaBeamformingImt(par,300,-10)
        par = self.param.get_antenna_parameters("UE","TX")
        self.antenna2 = AntennaBeamformingImt(par,-33.21,-5.31)

    def test_azimuth(self):
        self.assertEqual(self.antenna1.azimuth,300)
        self.assertEqual(self.antenna2.azimuth,-33.21)

    def test_elevation(self):
        self.assertEqual(self.antenna1.elevation,-10)
        self.assertEqual(self.antenna2.elevation,-5.31)

    def test_g_max(self):
        self.assertEqual(self.antenna1.element.g_max,5)
        self.assertEqual(self.antenna2.element.g_max,10)

    def test_phi_3db(self):
        self.assertEqual(self.antenna1.element.phi_deg_3db,80)
        self.assertEqual(self.antenna2.element.phi_deg_3db,75)

    def test_theta_3db(self):
        self.assertEqual(self.antenna1.element.theta_deg_3db,60)
        self.assertEqual(self.antenna2.element.theta_deg_3db,65)

    def test_am(self):
        self.assertEqual(self.antenna1.element.am,30)
        self.assertEqual(self.antenna2.element.am,25)

    def test_sla_v(self):
        self.assertEqual(self.antenna1.element.sla_v,30)
        self.assertEqual(self.antenna2.element.sla_v,35)

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

    def test_beams_list(self):
        self.assertEqual(len(self.antenna1.beams_list),0)
        self.assertEqual(len(self.antenna2.beams_list),0)

    def test_w_vec_list(self):
        self.assertEqual(len(self.antenna1.w_vec_list),0)
        self.assertEqual(len(self.antenna2.w_vec_list),0)

    def test_super_position_vector(self):
        # Error margin
        eps = 1e-5

        # Test 1
        phi = 0
        theta = 0
        v_vec = self.antenna2._super_position_vector(phi, theta)
        expected_v_vec = np.array([[1.0, 1.0],[-1.0, -1.0]])
        self.assertTrue(np.allclose(np.real(v_vec),\
                                    np.real(expected_v_vec),rtol = eps))
        self.assertTrue(np.allclose(np.imag(v_vec),\
                                    np.imag(expected_v_vec),rtol = eps))

        # Test 2
        phi = 90
        theta = 90
        v_vec = self.antenna2._super_position_vector(phi, theta)
        expected_v_vec = np.array([[1.0, -1.0],[1.0, -1.0]])
        self.assertTrue(np.allclose(np.real(v_vec),\
                                    np.real(expected_v_vec),rtol = eps))
        self.assertTrue(np.allclose(np.imag(v_vec),\
                                    np.imag(expected_v_vec),rtol = eps))

        # Test 3
        phi = 45
        theta = 45
        v_vec = self.antenna2._super_position_vector(phi, theta)
        expected_v_vec = np.array([[1.0 + 0.0j, 0.0 + 1.0j],\
                    [-0.6056998+0.7956932j, -0.7956932-0.6056998j]])
        self.assertTrue(np.allclose(np.real(v_vec),\
                                    np.real(expected_v_vec),rtol = eps))
        self.assertTrue(np.allclose(np.imag(v_vec),\
                                    np.imag(expected_v_vec),rtol = eps))

        # Test 4
        phi = 60
        theta = 90
        v_vec = self.antenna2._super_position_vector(phi, theta)
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
        w_vec = self.antenna2._weight_vector(phi_scan, theta_tilt)
        expected_w_vec = np.array([[0.5, 0.5],[0.5, 0.5]])
        self.assertTrue(np.allclose(np.real(w_vec),\
                                    np.real(expected_w_vec),rtol = eps))
        self.assertTrue(np.allclose(np.imag(w_vec),\
                                    np.imag(expected_w_vec),rtol = eps))

        # Test 2
        phi_scan = 90
        theta_tilt = 90
        w_vec = self.antenna2._weight_vector(phi_scan, theta_tilt)
        expected_w_vec = np.array([[0.5, 0.5],[-0.5, -0.5]])
        self.assertTrue(np.allclose(np.real(w_vec),\
                                    np.real(expected_w_vec),rtol = eps))
        self.assertTrue(np.allclose(np.imag(w_vec),\
                                    np.imag(expected_w_vec),rtol = eps))

        # Test 3
        phi_scan = 45
        theta_tilt = 45
        w_vec = self.antenna2._weight_vector(phi_scan, theta_tilt)
        expected_w_vec = np.array([[0.5 + 0.0j, 0.0 - 0.5j],\
                    [-0.3028499+0.3978466j, 0.3978466+0.3028499j]])
        self.assertTrue(np.allclose(np.real(w_vec),\
                                    np.real(expected_w_vec),rtol = eps))
        self.assertTrue(np.allclose(np.imag(w_vec),\
                                    np.imag(expected_w_vec),rtol = eps))

        # Test 4
        phi_scan = 0
        theta_tilt = 90
        w_vec = self.antenna2._weight_vector(phi_scan, theta_tilt)
        expected_w_vec = np.array([[0.5, 0.5],[-0.5, -0.5]])
        self.assertTrue(np.allclose(np.real(w_vec),\
                                    np.real(expected_w_vec),rtol = eps))
        self.assertTrue(np.allclose(np.imag(w_vec),\
                                    np.imag(expected_w_vec),rtol = eps))

        # Test 5
        phi_scan = 45
        theta_tilt = 30
        w_vec = self.antenna2._weight_vector(phi_scan, theta_tilt)
        expected_w_vec = np.array([[0.5 + 0.0j, -0.172870 - 0.469169j],\
                                   [0.0 + 0.5j,  0.469165 - 0.172870j]])
        self.assertTrue(np.allclose(np.real(w_vec),\
                                    np.real(expected_w_vec),rtol = eps))
        self.assertTrue(np.allclose(np.imag(w_vec),\
                                    np.imag(expected_w_vec),rtol = eps))


    def test_add_beam(self):
        # Error margin and antenna object
        eps = 1e-5
        par = self.param.get_antenna_parameters("UE","TX")
        self.antenna2 = AntennaBeamformingImt(par,0,0)

        # Add first beam
        phi_scan = 45
        theta_tilt = 120
        self.antenna2.add_beam(phi_scan,theta_tilt)

        self.assertEqual(len(self.antenna2.beams_list),1)
        self.assertEqual(len(self.antenna2.w_vec_list),1)

        # Add second beam
        phi_scan = 90
        theta_tilt = 180.0
        self.antenna2.add_beam(phi_scan,theta_tilt)

        # Test beams_list
        self.assertEqual(len(self.antenna2.beams_list),2)
        self.assertEqual(len(self.antenna2.w_vec_list),2)

        # Test first beam
        self.assertAlmostEqual(self.antenna2.beams_list[0][0],45,delta=eps)
        self.assertAlmostEqual(self.antenna2.beams_list[0][1],30,delta=eps)

        w_vec = self.antenna2.w_vec_list[0]
        expected_w_vec = np.array([[0.5 + 0.0j, -0.172870 - 0.469169j],\
                                   [0.0 + 0.5j,  0.469165 - 0.172870j]])
        self.assertTrue(np.allclose(np.real(w_vec),\
                                    np.real(expected_w_vec),rtol = eps))
        self.assertTrue(np.allclose(np.imag(w_vec),\
                                    np.imag(expected_w_vec),rtol = eps))

        # Test second beam
        self.assertEqual(self.antenna2.beams_list[1][0],90)
        self.assertEqual(self.antenna2.beams_list[1][1],90)

        w_vec = self.antenna2.w_vec_list[1]
        expected_w_vec = np.array([[0.5, 0.5],[-0.5, -0.5]])
        self.assertTrue(np.allclose(np.real(w_vec),\
                                    np.real(expected_w_vec),rtol = eps))
        self.assertTrue(np.allclose(np.imag(w_vec),\
                                    np.imag(expected_w_vec),rtol = eps))

        # Reset beams and test
        self.antenna2.reset_beams()
        self.assertEqual(len(self.antenna2.beams_list),0)
        self.assertEqual(len(self.antenna2.w_vec_list),0)

    def test_beam_gain(self):
        # Error margin and antenna
        eps = 1e-4
        par = self.param.get_antenna_parameters("UE","TX")
        self.antenna2 = AntennaBeamformingImt(par,0,0)

        # Test 1
        phi = 45
        theta = 45
        beam = 0
        phi_scan = 45
        theta_tilt = 135
        self.antenna2.add_beam(phi_scan,theta_tilt)
        beam_g = self.antenna2._beam_gain(phi,theta,beam)
        self.assertAlmostEqual(beam_g,1.594268,delta = eps)

        # Test 2
        phi = 0
        theta = 60
        beam = 1
        phi_scan = 45
        theta_tilt = 180
        self.antenna2.add_beam(phi_scan,theta_tilt)
        beam_g = self.antenna2._beam_gain(phi,theta,beam)
        self.assertAlmostEqual(beam_g,10.454087,delta = eps)

        # Test 3
        phi = 32.5
        theta = 115.2
        beam_g = self.antenna2._beam_gain(phi,theta)
        self.assertAlmostEqual(beam_g,11.9636,delta = eps)

    def test_calculate_gain(self):
        # Error margin and antenna
        eps = 1e-4
        par = self.param.get_antenna_parameters("BS","RX")
        self.antenna1 = AntennaBeamformingImt(par,0,0)
        par = self.param.get_antenna_parameters("UE","TX")
        self.antenna2 = AntennaBeamformingImt(par,0,0)

        # Test 1
        phi_vec = np.array([45.0, 32.5])
        theta_vec = np.array([45.0, 115.2])
        gains = self.antenna2.calculate_gain(phi_vec=phi_vec, theta_vec=theta_vec)
        npt.assert_allclose(gains,np.array([5.9491,11.9636]),atol=eps)

        # Test 2
        phi = 0.0
        theta = 60.0
        phi_scan = 45
        theta_tilt = 180
        self.antenna2.add_beam(phi_scan, theta_tilt)
        beams_l = np.zeros_like(phi_vec, dtype=int)
        gains = self.antenna2.calculate_gain(phi_vec=phi, theta_vec=theta,
                                             beams_l=beams_l)
        npt.assert_allclose(gains,np.array([10.454087]),atol=eps)

        # Test 3
        phi = 40
        theta = 100
        gains = self.antenna1.calculate_gain(phi_vec=phi, theta_vec=theta,
                                             co_channel=False)
        npt.assert_allclose(gains,np.array([1.6667]),atol=eps)
        
    def test_normalization(self):
        # Create dummy normalization data
        normalization = True
        norm_data = {'norm_file': 'dummy_file.npz',
                     'resolution': 1,
                     'phi_range': (-180,+180),
                     'theta_range': (0,180),
                     'correction_factor_co_channel': np.ones((360,180)),
                     'error_co_channel': 0.0,
                     'correction_factor_adj_channel': 5,
                     'error_adj_channel': 0.0,
                     'parameters': None}
        element_pattern = "M2101"
        element_max_g = 5
        element_phi_deg_3db = 65
        element_theta_deg_3db = 65
        element_am = 30
        element_sla_v = 30
        n_rows = 8
        n_columns = 8
        horiz_spacing = 0.5
        vert_spacing = 0.5
        down_tilt = 0
        par = AntennaPar(normalization,
                         norm_data,
                         element_pattern,
                         element_max_g,
                         element_phi_deg_3db,
                         element_theta_deg_3db,
                         element_am,
                         element_sla_v,
                         n_rows,
                         n_columns,
                         horiz_spacing,
                         vert_spacing,
                         down_tilt)
        
        
        # Create antenna objects
        self.antenna3 = AntennaBeamformingImt(par,0.0,0.0) # Normalized
        par = par._replace(normalization = False)
        self.antenna4 = AntennaBeamformingImt(par,0.0,0.0) # Unormalized
        
        # Test co-channel gain: no beam
        phi_v = np.array([11.79, -0.71])
        theta_v = np.array([50.31, 120.51])
        gain_ref = self.antenna4.calculate_gain(phi_vec=phi_v, theta_vec=theta_v)
        gain = self.antenna3.calculate_gain(phi_vec=phi_v, theta_vec=theta_v)
        npt.assert_equal(gain,gain_ref + 1)
        
        # Test co-channel gain: add beam
        phi_scan = 11.79
        theta_tilt = 185.31
        self.antenna3.add_beam(phi_scan, theta_tilt)
        self.antenna4.add_beam(phi_scan, theta_tilt)
        beams_l = np.zeros_like(phi_v, dtype=int)
        gain_ref = self.antenna4.calculate_gain(phi_vec=phi_v, theta_vec=theta_v,
                                             beams_l=beams_l)
        gain = self.antenna3.calculate_gain(phi_vec=phi_v, theta_vec=theta_v,
                                             beams_l=beams_l)
        npt.assert_equal(gain,gain_ref + 1)
        
        # Test adjacent channel
        phi_v = np.array([11.79, -0.71])
        theta_v = np.array([50.31, 120.51])
        gain_ref = self.antenna4.calculate_gain(phi_vec=phi_v, theta_vec=theta_v, co_channel=False)
        gain = self.antenna3.calculate_gain(phi_vec=phi_v, theta_vec=theta_v, co_channel=False)
        npt.assert_equal(gain,gain_ref + 5)

    def test_to_local_coord(self):
        # Test 1
        # Create antenna object
        par = self.param.get_antenna_parameters("BS","RX")
        azi = 0.0
        ele = 90.0
        self.antenna3 = AntennaBeamformingImt(par,azi,ele)
        
        # Angles to be converted: z, x and y axis
        phi = np.array([   0, 90, 0])
        theta = np.array([90, 90, 0])

        # Convert to local coordinates
        lo_phi, lo_theta = self.antenna3.to_local_coord(phi, theta)
        exp_lo_phi = np.array([0,90,0])
        exp_lo_theta = np.array([180,90,90])
        npt.assert_array_almost_equal(lo_phi,exp_lo_phi,decimal=2)
        npt.assert_array_almost_equal(lo_theta,exp_lo_theta,decimal=2)
        
        # Test 2
        # Create antenna object
        azi = -15.0
        ele = 90.0
        self.antenna3 = AntennaBeamformingImt(par,azi,ele)
        
        # Angles to be converted
        phi = np.array([  0,  0, 90, 10])
        theta = np.array([0, 90, 90, 90])

        # Convert to local coordinates
        lo_phi, lo_theta = self.antenna3.to_local_coord(phi, theta)
        exp_lo_phi = np.array([0, 90, 90, 90])
        exp_lo_theta = np.array([90, 165, 75, 155])
        npt.assert_array_almost_equal(lo_phi,exp_lo_phi,decimal=2)
        npt.assert_array_almost_equal(lo_theta,exp_lo_theta,decimal=2)
        
        # Test 3
        # Create antenna object
        azi = +15.0
        ele = 90.0
        self.antenna3 = AntennaBeamformingImt(par,azi,ele)
        
        # Angles to be converted: z, x and y axis
        phi = np.array([  0,  0, 90])
        theta = np.array([0, 90, 90])

        # Convert to local coordinates
        lo_phi, lo_theta = self.antenna3.to_local_coord(phi, theta)
        exp_lo_phi = np.array([0,-90, 90])
        exp_lo_theta = np.array([90, 165, 105])
        npt.assert_array_almost_equal(lo_phi,exp_lo_phi,decimal=2)
        npt.assert_array_almost_equal(lo_theta,exp_lo_theta,decimal=2)
        
        # Test 4
        # Create antenna object
        azi = 0.0
        ele = 15.0
        self.antenna3 = AntennaBeamformingImt(par,azi,ele)
        
        # Angles to be converted: z, x and y axis
        phi = np.array([  0,  0, 90])
        theta = np.array([0, 90, 90])

        # Convert to local coordinates
        lo_phi, lo_theta = self.antenna3.to_local_coord(phi, theta)
        exp_lo_phi = np.array([0, 0, 90])
        exp_lo_theta = np.array([15, 105, 90])
        npt.assert_array_almost_equal(lo_phi,exp_lo_phi,decimal=2)
        npt.assert_array_almost_equal(lo_theta,exp_lo_theta,decimal=2)
        
        # Test 5
        # Create antenna object
        azi = 0.0
        ele = -15.0
        self.antenna3 = AntennaBeamformingImt(par,azi,ele)
        
        # Angles to be converted: z, x and y axis
        phi = np.array([  0,  0, 90])
        theta = np.array([0, 90, 90])

        # Convert to local coordinates
        lo_phi, lo_theta = self.antenna3.to_local_coord(phi, theta)
        exp_lo_phi = np.array([180, 0, 90])
        exp_lo_theta = np.array([15, 75, 90])
        npt.assert_array_almost_equal(lo_phi,exp_lo_phi,decimal=2)
        npt.assert_array_almost_equal(lo_theta,exp_lo_theta,decimal=2)
        
        # Test 6
        # Create antenna object
        azi = 180.0
        ele = 45.0
        self.antenna3 = AntennaBeamformingImt(par,azi,ele)
        
        # Angles to be converted: z, x and y axis
        phi = np.array([  0,  0, 90])
        theta = np.array([0, 90, 90])

        # Convert to local coordinates
        lo_phi, lo_theta = self.antenna3.to_local_coord(phi, theta)
        exp_lo_phi = np.array([0, -180, -90])
        exp_lo_theta = np.array([45, 45, 90])
        npt.assert_array_almost_equal(lo_phi,exp_lo_phi,decimal=2)
        npt.assert_array_almost_equal(lo_theta,exp_lo_theta,decimal=2)
        
        
if __name__ == '__main__':
    unittest.main()
#    
#    suite = unittest.TestSuite()
#    suite.addTest(AntennaBeamformingImtTest('test_calculate_gain'))
#    unittest.TextTestRunner().run(suite)
