# -*- coding: utf-8 -*-
"""
Created on Mon Feb  6 18:13:11 2017

@author: edgar
"""

import unittest
import numpy as np
import numpy.testing as npt

from sharc.support.enumerations import StationType
from sharc.antenna.antenna_beamforming_imt import AntennaBeamformingImt
from sharc.parameters.parameters_antenna_imt import ParametersAntennaImt
from sharc.station import Station
from sharc.station_manager import StationManager


class StationManagerTest(unittest.TestCase):

    def setUp(self):
        #Array parameters
        self.param = ParametersAntennaImt()

        self.param.bs_element_pattern = "M2101"
        self.param.bs_downtilt_deg = 0

        self.param.bs_rx_element_max_g = 10
        self.param.bs_rx_element_phi_deg_3db = 65
        self.param.bs_rx_element_theta_deg_3db = 75
        self.param.bs_rx_element_am = 35
        self.param.bs_rx_element_sla_v = 25
        self.param.bs_rx_n_rows = 8
        self.param.bs_rx_n_columns = 8
        self.param.bs_rx_element_horiz_spacing = 0.5
        self.param.bs_rx_element_vert_spacing = 0.5

        self.param.bs_tx_element_max_g = 5
        self.param.bs_tx_element_phi_deg_3db = 80
        self.param.bs_tx_element_theta_deg_3db = 60
        self.param.bs_tx_element_am = 30
        self.param.bs_tx_element_sla_v = 30
        self.param.bs_tx_n_rows = 16
        self.param.bs_tx_n_columns = 16
        self.param.bs_tx_element_horiz_spacing = 1
        self.param.bs_tx_element_vert_spacing = 1

        self.param.ue_element_pattern = "M2101"

        self.param.ue_rx_element_max_g = 10
        self.param.ue_rx_element_phi_deg_3db = 75
        self.param.ue_rx_element_theta_deg_3db = 65
        self.param.ue_rx_element_am = 25
        self.param.ue_rx_element_sla_v = 35
        self.param.ue_rx_n_rows = 2
        self.param.ue_rx_n_columns = 2
        self.param.ue_rx_element_horiz_spacing = 0.5
        self.param.ue_rx_element_vert_spacing = 0.5

        self.param.ue_tx_element_max_g = 10
        self.param.ue_tx_element_phi_deg_3db = 75
        self.param.ue_tx_element_theta_deg_3db = 65
        self.param.ue_tx_element_am = 25
        self.param.ue_tx_element_sla_v = 35
        self.param.ue_tx_n_rows = 2
        self.param.ue_tx_n_columns = 2
        self.param.ue_tx_element_horiz_spacing = 0.5
        self.param.ue_tx_element_vert_spacing = 0.5

        self.station_manager = StationManager(3)
        self.station_manager.x = np.array([10, 20, 30])
        self.station_manager.y = np.array([15, 25, 35])
        self.station_manager.height = np.array([1, 2, 3])
        # this is for downlink
        self.station_manager.tx_power = dict({0: [27, 30], 1: [35], 2: [40]})
        self.station_manager.rx_power = np.array([-50, -35, -10])
        par = self.param.get_antenna_parameters("BS","TX")
        self.station_manager.antenna = np.array([AntennaBeamformingImt(par,60,-10), AntennaBeamformingImt(par,180,-10), AntennaBeamformingImt(par,300,-10)])
        self.station_manager.station_type = StationType.IMT_BS

        self.station_manager2 = StationManager(2)
        self.station_manager2.x = np.array([100, 200])
        self.station_manager2.y = np.array([105, 250])
        self.station_manager2.height = np.array([4, 5])
        # this is for downlink
        self.station_manager2.tx_power = dict({0: [25], 1: [28,35]})
        self.station_manager2.rx_power = np.array([-50, -35])
        par = self.param.get_antenna_parameters("BS","RX")
        self.station_manager2.antenna = np.array([AntennaBeamformingImt(par,0,-5), AntennaBeamformingImt(par,180,-5)])
        self.station_manager2.station_type = StationType.IMT_BS

        self.station_manager3 = StationManager(1)
        self.station_manager3.x = np.array([300])
        self.station_manager3.y = np.array([400])
        self.station_manager3.height = np.array([2])
        # this is for uplink
        self.station_manager3.tx_power = 22
        self.station_manager3.rx_power = np.array([-50,-35])
        par = self.param.get_antenna_parameters("UE","TX")
        self.station_manager3.antenna = np.array([AntennaBeamformingImt(par,0,-30), AntennaBeamformingImt(par,35,45)])
        self.station_manager3.station_type = StationType.IMT_UE

        self.station = Station()
        self.station.id = 0
        self.station.x = 10
        self.station.y = 15
        self.station.height = 1
        self.station.tx_power = 30
        self.station.rx_power = -50
        par = self.param.get_antenna_parameters("UE","TX")
        self.station.antenna = AntennaBeamformingImt(par,100,-10)
        self.station.station_type = StationType.IMT_UE

        self.station2 = Station()
        self.station2.id = 1
        self.station2.x = 20
        self.station2.y = 25
        self.station2.height = 2
        self.station2.tx_power = 35
        self.station2.rx_power = -35
        par = self.param.get_antenna_parameters("BS","RX")
        self.station2.antenna = AntennaBeamformingImt(par,-90,-15)
        self.station2.station_type = StationType.IMT_BS


    def test_num_stations(self):
        self.assertEqual(self.station_manager.num_stations, 3)
        self.assertEqual(self.station_manager2.num_stations, 2)
        self.assertEqual(self.station_manager3.num_stations, 1)

    def test_station_type(self):
        self.assertEqual(self.station_manager.station_type, StationType.IMT_BS)
        self.assertEqual(self.station_manager2.station_type, StationType.IMT_BS)
        self.assertEqual(self.station_manager3.station_type, StationType.IMT_UE)

    def test_x(self):
        # get a single value from the original array
        self.assertEqual(self.station_manager.x[0], 10)
        # get two specific values
        npt.assert_array_equal(self.station_manager.x[[1,2]], [20,30])
        # get values in reverse order
        npt.assert_array_equal(self.station_manager.x[[2,1,0]], [30,20,10])
        # get all values (no need to specify the id's)
        npt.assert_array_equal(self.station_manager.x, [10,20,30])
        # set a single value and get it
        self.station_manager.x[0] = 8
        npt.assert_array_equal(self.station_manager.x[[0,1]], [8,20])
        # set two values and then get all values
        self.station_manager.x[[1,2]] = [16,32]
        npt.assert_array_equal(self.station_manager.x, [8,16,32])

    def test_y(self):
        # get a single value from the original array
        self.assertEqual(self.station_manager.y[0], 15)
        # get two specific values
        npt.assert_array_equal(self.station_manager.y[[1,2]], [25,35])
        # get values in reverse order
        npt.assert_array_equal(self.station_manager.y[[2,1,0]], [35,25,15])
        # get all values (no need to specify the id's)
        npt.assert_array_equal(self.station_manager.y, [15,25,35])
        # set a single value and get it
        self.station_manager.y[1] = 9
        npt.assert_array_equal(self.station_manager.y[[0,1]], [15,9])
        # set two values and then get all values
        self.station_manager.y[[0,2]] = [7,21]
        npt.assert_array_equal(self.station_manager.y, [7,9,21])

    def test_height(self):
        # get a single value from the original array
        self.assertEqual(self.station_manager.height[0], 1)
        # get two specific values
        npt.assert_array_equal(self.station_manager.height[[0,2]], [1,3])
        # get values in reverse order
        npt.assert_array_equal(self.station_manager.height[[2,1,0]], [3,2,1])
        # get all values (no need to specify the id's)
        npt.assert_array_equal(self.station_manager.height, [1,2,3])
        # set a single value and get it
        self.station_manager.height[1] = 7
        npt.assert_array_equal(self.station_manager.height[[1,2]], [7,3])
        # set two values and then get all values
        self.station_manager.height[[0,2]] = [5,4]
        npt.assert_array_equal(self.station_manager.height, [5,7,4])

    def test_tx_power(self):
        # get a single value from the original array
        self.assertEqual(self.station_manager.tx_power[0], [27,30])
        self.assertEqual(self.station_manager.tx_power[1], [35])
        # get all values (no need to specify the id's)
        npt.assert_array_equal(self.station_manager.tx_power, dict({0: [27, 30], 1: [35], 2: [40]}))
        # set a single value and get it
        self.station_manager.tx_power[0] = [33,38]
        npt.assert_array_equal(self.station_manager.tx_power[0], [33,38])
        # set two values and then get all values
        self.station_manager.tx_power[2] = [20,25]
        npt.assert_array_equal(self.station_manager.tx_power, dict({0: [33,38], 1: [35], 2: [20,25]}))

    def test_rx_power(self):
        # get a single value from the original array
        self.assertEqual(self.station_manager.rx_power[2], -10)
        # get two specific values
        npt.assert_array_equal(self.station_manager.rx_power[[0,1]], [-50,-35])
        # get values in reverse order
        npt.assert_array_equal(self.station_manager.rx_power[[2,1,0]], [-10,-35,-50] )
        # get all values (no need to specify the id's)
        npt.assert_array_equal(self.station_manager.rx_power, [-50,-35,-10])
        # set a single value and get it
        self.station_manager.rx_power[2] = -15
        npt.assert_array_equal(self.station_manager.rx_power[[2,0]], [-15,-50])
        # set two values and then get all values
        self.station_manager.rx_power[[0,1]] = [-60,-30]
        npt.assert_array_equal(self.station_manager.rx_power, [-60,-30,-15])

    def test_antenna(self):
        self.assertEqual(self.station_manager.antenna[0].azimuth, 60)
        self.assertEqual(self.station_manager.antenna[0].elevation, -10)
        self.assertEqual(self.station_manager.antenna[0].n_rows, 16)
        self.assertEqual(self.station_manager.antenna[0].n_cols, 16)
        self.assertEqual(self.station_manager.antenna[0].dh, 1)
        self.assertEqual(self.station_manager.antenna[0].dv, 1)
        self.assertEqual(self.station_manager.antenna[0].element.g_max, 5)
        self.assertEqual(self.station_manager.antenna[0].element.phi_deg_3db, 80)
        self.assertEqual(self.station_manager.antenna[0].element.theta_deg_3db, 60)
        self.assertEqual(self.station_manager.antenna[0].element.am, 30)
        self.assertEqual(self.station_manager.antenna[0].element.sla_v, 30)

        self.assertEqual(self.station_manager.antenna[1].azimuth, 180)
        self.assertEqual(self.station_manager.antenna[1].elevation, -10)
        self.assertEqual(self.station_manager.antenna[1].n_rows, 16)
        self.assertEqual(self.station_manager.antenna[1].n_cols, 16)
        self.assertEqual(self.station_manager.antenna[1].dh, 1)
        self.assertEqual(self.station_manager.antenna[1].dv, 1)
        self.assertEqual(self.station_manager.antenna[1].element.g_max, 5)
        self.assertEqual(self.station_manager.antenna[1].element.phi_deg_3db, 80)
        self.assertEqual(self.station_manager.antenna[1].element.theta_deg_3db, 60)
        self.assertEqual(self.station_manager.antenna[1].element.am, 30)
        self.assertEqual(self.station_manager.antenna[1].element.sla_v, 30)

        self.assertEqual(self.station_manager.antenna[2].azimuth, 300)
        self.assertEqual(self.station_manager.antenna[2].elevation, -10)
        self.assertEqual(self.station_manager.antenna[2].n_rows, 16)
        self.assertEqual(self.station_manager.antenna[2].n_cols, 16)
        self.assertEqual(self.station_manager.antenna[2].dh, 1)
        self.assertEqual(self.station_manager.antenna[2].dv, 1)
        self.assertEqual(self.station_manager.antenna[2].element.g_max, 5)
        self.assertEqual(self.station_manager.antenna[2].element.phi_deg_3db, 80)
        self.assertEqual(self.station_manager.antenna[2].element.theta_deg_3db, 60)
        self.assertEqual(self.station_manager.antenna[2].element.am, 30)
        self.assertEqual(self.station_manager.antenna[2].element.sla_v, 30)

        par = self.param.get_antenna_parameters("BS","TX")
        self.station_manager.antenna[[0,2]] = np.array([AntennaBeamformingImt(par,0,-5), AntennaBeamformingImt(par,180,-5)])

        self.assertEqual(self.station_manager.antenna[0].azimuth, 0)
        self.assertEqual(self.station_manager.antenna[0].elevation, -5)
        self.assertEqual(self.station_manager.antenna[0].n_rows, 16)
        self.assertEqual(self.station_manager.antenna[0].n_cols, 16)
        self.assertEqual(self.station_manager.antenna[0].dh, 1)
        self.assertEqual(self.station_manager.antenna[0].dv, 1)
        self.assertEqual(self.station_manager.antenna[0].element.g_max, 5)
        self.assertEqual(self.station_manager.antenna[0].element.phi_deg_3db, 80)
        self.assertEqual(self.station_manager.antenna[0].element.theta_deg_3db, 60)
        self.assertEqual(self.station_manager.antenna[0].element.am, 30)
        self.assertEqual(self.station_manager.antenna[0].element.sla_v, 30)

        self.assertEqual(self.station_manager.antenna[2].azimuth, 180)
        self.assertEqual(self.station_manager.antenna[2].elevation, -5)
        self.assertEqual(self.station_manager.antenna[2].n_rows, 16)
        self.assertEqual(self.station_manager.antenna[2].n_cols, 16)
        self.assertEqual(self.station_manager.antenna[2].dh, 1)
        self.assertEqual(self.station_manager.antenna[2].dv, 1)
        self.assertEqual(self.station_manager.antenna[2].element.g_max, 5)
        self.assertEqual(self.station_manager.antenna[2].element.phi_deg_3db, 80)
        self.assertEqual(self.station_manager.antenna[2].element.theta_deg_3db, 60)
        self.assertEqual(self.station_manager.antenna[2].element.am, 30)
        self.assertEqual(self.station_manager.antenna[2].element.sla_v, 30)

    def test_station(self):
        # test if manager returns the correct station
        self.assertTrue(self.station == self.station_manager.get_station(0))
        self.assertTrue(self.station2 == self.station_manager.get_station(1))
        # now we change station id and verify if stations became different
        self.station.id = 1
        self.assertTrue(self.station != self.station_manager.get_station(0))
        # test station type
        self.assertEqual(self.station_manager.get_station(0).station_type,\
                         StationType.IMT_BS)

    def test_station_list(self):
        # test if manager returns the correct station list
        station_list = self.station_manager.get_station_list()
        self.assertTrue(self.station in station_list)
        self.assertTrue(self.station2 in station_list)
        #
        station_list = self.station_manager.get_station_list([0,2])
        self.assertTrue(self.station in station_list)
        self.assertTrue(self.station2 not in station_list)
        #
        station_list = self.station_manager.get_station_list([2])
        self.assertTrue(self.station not in station_list)
        self.assertTrue(self.station2 not in station_list)

    def test_distance_to(self):
        ref_distance = np.array([[ 356.405,  180.277]])
        distance = self.station_manager3.get_distance_to(self.station_manager2)
        npt.assert_allclose(distance, ref_distance, atol=1e-2)

        ref_distance = np.asarray([[ 127.279,  302.200],
                                   [ 113.137,  288.140],
                                   [  98.994,  274.089]])
        distance = self.station_manager.get_distance_to(self.station_manager2)
        npt.assert_allclose(distance, ref_distance, atol=1e-2)

    def test_3d_distance_to(self):
        ref_distance = np.asarray([[ 356.411,  180.302]])
        distance = self.station_manager3.get_3d_distance_to(self.station_manager2)
        npt.assert_allclose(distance, ref_distance, atol=1e-2)

        ref_distance = np.asarray([[ 127.314,  302.226],
                                   [ 113.154,  288.156],
                                   [  99,  274.096]])
        distance = self.station_manager.get_3d_distance_to(self.station_manager2)
        npt.assert_allclose(distance, ref_distance, atol=1e-2)

    def test_pointing_vector_to(self):
        eps = 1e-1
        # Test 1
        phi, theta = self.station_manager.get_pointing_vector_to(self.station_manager2)
        npt.assert_allclose(phi,np.array([[45.00, 51.04],
                                          [45.00, 51.34],
                                          [45.00, 51.67]]),atol=eps)
        npt.assert_allclose(theta,np.array([[88.65, 89.24],
                                            [88.89, 89.40],
                                            [89.42, 89.58]]),atol=eps)

        # Test 2
        phi, theta = self.station_manager2.get_pointing_vector_to(self.station_manager)
        npt.assert_allclose(phi,np.array([[-135.00, -135.00, -135.00],
                                          [-128.96, -128.66, -128.33]]),atol=eps)
        npt.assert_allclose(theta,np.array([[91.35, 91.01, 90.58],
                                            [90.76, 90.60, 90.42]]),atol=eps)

        # Test 3
        phi, theta = self.station_manager3.get_pointing_vector_to(self.station_manager2)
        npt.assert_allclose(phi,np.array([[-124.13, -123.69]]),atol=eps)
        npt.assert_allclose(theta,np.array([[89.73, 89.05]]),atol=eps)

        # Test 4
        phi, theta = self.station_manager2.get_pointing_vector_to(self.station_manager3)
        npt.assert_allclose(phi,np.array([[55.86], [56.31]]),atol=eps)
        npt.assert_allclose(theta,np.array([[90.32], [90.95]]),atol=eps)

    def test_off_axis_angle(self):
        sm1 = StationManager(1)
        sm1.x = np.array([0])
        sm1.y = np.array([0])
        sm1.height = np.array([0])
        sm1.azimuth = np.array([0])
        sm1.elevation = np.array([0])

        sm2 = StationManager(6)
        sm2.x = np.array([100, 100, 0, 100, 100, 100])
        sm2.y = np.array([0, 0, 100, 100, 100, 100])
        sm2.height = np.array([0, 100, 0, 0, 100, 100])
        sm2.azimuth = np.array([180, 180, 180, 180, 180, 225])
        sm2.elevation = np.array([0, 0, 0, 0, 0, 0])

        phi_ref = np.array([[0, 45, 90, 45, 54.73, 54.73]])
        npt.assert_allclose(phi_ref, sm1.get_off_axis_angle(sm2), atol=1e-2)

        #######################################################################
        sm3 = StationManager(1)
        sm3.x = np.array([0])
        sm3.y = np.array([0])
        sm3.height = np.array([0])
        sm3.azimuth = np.array([45])
        sm3.elevation = np.array([0])

        sm4 = StationManager(2)
        sm4.x = np.array([100, 60])
        sm4.y = np.array([100, 80])
        sm4.height = np.array([100, 100])
        sm4.azimuth = np.array([180, 180])
        sm4.elevation = np.array([0, 0])

        phi_ref = np.array([[35.26, 45.57]])
        npt.assert_allclose(phi_ref, sm3.get_off_axis_angle(sm4), atol=1e-2)

        #######################################################################
        sm5 = StationManager(1)
        sm5.x = np.array([0])
        sm5.y = np.array([0])
        sm5.height = np.array([0])
        sm5.azimuth = np.array([0])
        sm5.elevation = np.array([45])

        sm6 = StationManager(2)
        sm6.x = np.array([100, 100])
        sm6.y = np.array([0, 0])
        sm6.height = np.array([100, 100])
        sm6.azimuth = np.array([180, 180])
        sm6.elevation = np.array([0, 0])

        phi_ref = np.array([[0, 0]])
        npt.assert_allclose(phi_ref, sm5.get_off_axis_angle(sm6), atol=1e-2)

        #######################################################################
        sm6 = StationManager(1)
        sm6.x = np.array([0])
        sm6.y = np.array([0])
        sm6.height = np.array([100])
        sm6.azimuth = np.array([0])
        sm6.elevation = np.array([270])

        sm7 = StationManager(2)
        sm7.x = np.array([0, 100])
        sm7.y = np.array([0, 0])
        sm7.height = np.array([0, 0])
        sm7.azimuth = np.array([180, 180])
        sm7.elevation = np.array([0, 0])

        phi_ref = np.array([[0, 45]])
        npt.assert_allclose(phi_ref, sm6.get_off_axis_angle(sm7), atol=1e-2)          
        
        
    def test_elevation(self):
        sm1 = StationManager(1)
        sm1.x = np.array([0])
        sm1.y = np.array([0])
        sm1.height = np.array([10])
        
        sm2 = StationManager(6)
        sm2.x =      np.array([10, 10, 0,  0, 30, 20])
        sm2.y =      np.array([ 0,  0, 5, 10, 30, 20])
        sm2.height = np.array([10, 20, 5,  0, 20, 20])
        
        elevation_ref = np.array([[0, 45, -45, -45, 13.26, 19.47]])
        npt.assert_allclose(elevation_ref, sm1.get_elevation(sm2), atol=1e-2)    
        
        #######################################################################
        sm3 = StationManager(2)
        sm3.x = np.array([0, 30])
        sm3.y = np.array([0, 0])
        sm3.height = np.array([10, 10])
        
        sm4 = StationManager(2)
        sm4.x =      np.array([10, 10])
        sm4.y =      np.array([ 0,  0])
        sm4.height = np.array([10, 20])
        
        elevation_ref = np.array([[0, 45], [0, 26.56]])
        npt.assert_allclose(elevation_ref, sm3.get_elevation(sm4), atol=1e-2)          

if __name__ == '__main__':
    unittest.main()
