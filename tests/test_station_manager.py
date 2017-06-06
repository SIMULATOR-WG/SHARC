# -*- coding: utf-8 -*-
"""
Created on Mon Feb  6 18:13:11 2017

@author: edgar
"""

import unittest
import numpy as np
import numpy.testing as npt

from sharc.antenna.antenna_omni import AntennaOmni
from sharc.station import Station
from sharc.station_manager import StationManager


class StationManagerTest(unittest.TestCase):
    
    def setUp(self):
        self.station_manager = StationManager(3)
        self.station_manager.x = [10, 20, 30]
        self.station_manager.y = [15, 25, 35]
        self.station_manager.height = [1, 2, 3]
        # this is for downlink
        self.station_manager.tx_power = dict({0: [27, 30], 1: [35], 2: [40]})
        self.station_manager.rx_power = [-50, -35, -10]
        self.station_manager.tx_antenna = [AntennaOmni(10), AntennaOmni(25), AntennaOmni(30)]
        self.station_manager.rx_antenna = [AntennaOmni(5), AntennaOmni(15), AntennaOmni(20)]
        
        self.station_manager2 = StationManager(2)
        self.station_manager2.x = [100, 200]
        self.station_manager2.y = [105, 250]
        self.station_manager2.height = [4, 5]
        # this is for downlink
        self.station_manager2.tx_power = dict({0: [25], 1: [28,35]})
        self.station_manager2.rx_power = [-50, -35]
        self.station_manager2.tx_antenna = [AntennaOmni(10), AntennaOmni(25)]
        self.station_manager2.rx_antenna = [AntennaOmni(5), AntennaOmni(15)]      
        
        self.station_manager3 = StationManager(1)
        self.station_manager3.x = [300]
        self.station_manager3.y = [400]
        self.station_manager3.height = [2]
        # this is for uplink
        self.station_manager3.tx_power = 22
        self.station_manager3.rx_power = [-50,-35]
        self.station_manager3.tx_antenna = [AntennaOmni(10), AntennaOmni(25)]
        self.station_manager3.rx_antenna = [AntennaOmni(5), AntennaOmni(15)]
        
        self.station = Station()
        self.station.id = 0
        self.station.x = 10
        self.station.y = 15
        self.station.height = 1
        self.station.tx_power = 30
        self.station.rx_power = -50
        self.station.tx_antenna = AntennaOmni(10)
        self.station.rx_antenna = AntennaOmni(5)

        self.station2 = Station()
        self.station2.id = 1
        self.station2.x = 20
        self.station2.y = 25
        self.station2.height = 2
        self.station2.tx_power = 35
        self.station2.rx_power = -35
        self.station2.tx_antenna = AntennaOmni(25)
        self.station2.rx_antenna = AntennaOmni(15)
        
        
    def test_num_stations(self):
        self.assertEqual(self.station_manager.num_stations, 3)
        self.assertEqual(self.station_manager2.num_stations, 2)
        self.assertEqual(self.station_manager3.num_stations, 1)

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
        
    def test_tx_antenna(self):
        self.assertEqual(self.station_manager.tx_antenna[0], 10)
        self.assertEqual(self.station_manager.tx_antenna[1], 25)
        self.assertEqual(self.station_manager.tx_antenna[2], AntennaOmni(30))
        npt.assert_array_equal(self.station_manager.tx_antenna[[0,1]], [10,25])
        npt.assert_array_equal(self.station_manager.tx_antenna[[0,2]], [10,30])
        npt.assert_array_equal(self.station_manager.tx_antenna[[2,1]], [30,25])
        npt.assert_array_equal(self.station_manager.tx_antenna[[2,1,0]], [30,25,10])
        npt.assert_array_equal(self.station_manager.tx_antenna, [10,25,30])
        npt.assert_array_equal(self.station_manager.tx_antenna, [AntennaOmni(10),AntennaOmni(25),30])
        # test setting one antenna
        self.station_manager.tx_antenna[0] = AntennaOmni(2)
        npt.assert_array_equal(self.station_manager.tx_antenna, [2,25,30])
        # test setting two antennas
        self.station_manager.tx_antenna[[1,0]] = [AntennaOmni(2),AntennaOmni(4)]
        npt.assert_array_equal(self.station_manager.tx_antenna, [4,2,30])
        
    def test_rx_antenna(self):
        self.assertEqual(self.station_manager.rx_antenna[0], 5)
        self.assertEqual(self.station_manager.rx_antenna[1], 15)
        self.assertEqual(self.station_manager.rx_antenna[2], AntennaOmni(20))
        npt.assert_array_equal(self.station_manager.rx_antenna[[0,1]], [5,15])
        npt.assert_array_equal(self.station_manager.rx_antenna[[0,2]], [5,20])
        npt.assert_array_equal(self.station_manager.rx_antenna[[2,1]], [20,15])
        npt.assert_array_equal(self.station_manager.rx_antenna[[2,1,0]], [20,15,5])
        npt.assert_array_equal(self.station_manager.rx_antenna, [5,15,20])
        npt.assert_array_equal(self.station_manager.rx_antenna, [AntennaOmni(5),AntennaOmni(15),20])
        # test setting one antenna
        self.station_manager.rx_antenna[0] = AntennaOmni(7)
        npt.assert_array_equal(self.station_manager.rx_antenna, [7,15,20])
        # test setting two antennas
        self.station_manager.rx_antenna[[1,0]] = [AntennaOmni(8),AntennaOmni(9)]
        npt.assert_array_equal(self.station_manager.rx_antenna, [9,8,20])
        
    def test_station(self):
        # test if manager returns the correct station
        self.assertTrue(self.station == self.station_manager.get_station(0))
        self.assertTrue(self.station2 == self.station_manager.get_station(1))
        # now we change station id and verify if stations became different
        self.station.id = 1
        self.assertTrue(self.station != self.station_manager.get_station(0))

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
        
if __name__ == '__main__':
    unittest.main()
                