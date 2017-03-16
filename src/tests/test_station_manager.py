# -*- coding: utf-8 -*-
"""
Created on Mon Feb  6 18:13:11 2017

@author: edgar
"""

import unittest
import numpy as np

from antenna import Antenna
from station import Station
from station_manager import StationManager


class StationManagerTest(unittest.TestCase):
    """
    This test case contains sereval tests that could be written using the
    assert.allclose statement. But when we do it, it crashes Spyder. That is
    the reason for wrinting np.all(np.equal())
    """
    
    def setUp(self):
        self.station_manager = StationManager(3)
        self.station_manager.set_x([10,20,30])
        self.station_manager.set_y([15,25,35])
        self.station_manager.set_height([1,2,3])
        self.station_manager.set_tx_power([30,35,40])
        self.station_manager.set_rx_power([-50,-35,-10])
        self.station_manager.set_tx_antenna([Antenna(10),Antenna(25),Antenna(30)])
        self.station_manager.set_rx_antenna([Antenna(5),Antenna(15),Antenna(20)])
        
        self.station_manager2 = StationManager(2)
        self.station_manager2.set_x([100,200])
        self.station_manager2.set_y([105,250])
        self.station_manager2.set_height([4,5])
        self.station_manager2.set_tx_power([20,25])
        self.station_manager2.set_rx_power([-50,-35])
        self.station_manager2.set_tx_antenna([Antenna(10),Antenna(25)])
        self.station_manager2.set_rx_antenna([Antenna(5),Antenna(15)])        
        
        self.station_manager3 = StationManager(1)
        self.station_manager3.set_x([300])
        self.station_manager3.set_y([400])
        self.station_manager3.set_height([2])
        self.station_manager3.set_tx_power([22])
        self.station_manager3.set_rx_power([-50,-35])
        self.station_manager3.set_tx_antenna([Antenna(10),Antenna(25)])
        self.station_manager3.set_rx_antenna([Antenna(5),Antenna(15)])         
        
        self.station = Station()
        self.station.set_id(0)
        self.station.set_x(10)
        self.station.set_y(15)
        self.station.set_height(1)
        self.station.set_tx_power(30)
        self.station.set_rx_power(-50)
        self.station.set_tx_antenna(Antenna(10))
        self.station.set_rx_antenna(Antenna(5))        

        self.station2 = Station()
        self.station2.set_id(1)
        self.station2.set_x(20)
        self.station2.set_y(25)
        self.station2.set_height(2)
        self.station2.set_tx_power(35)
        self.station2.set_rx_power(-35)
        self.station2.set_tx_antenna(Antenna(25))
        self.station2.set_rx_antenna(Antenna(15))     
        
        
    def test_num_stations(self):
        self.assertEqual(3, self.station_manager.get_num_stations())
        self.station_manager.reset(5)
        self.assertEqual(5, self.station_manager.get_num_stations())

    def test_x(self):
        # get a single value from the original array
        self.assertEqual(10, self.station_manager.get_x(0))
        # get two specific values
        self.assertTrue(np.all(np.equal([20,30], self.station_manager.get_x([1,2]))))
        # get values in reverse order
        self.assertTrue(np.all(np.equal([30,20,10], self.station_manager.get_x([2,1,0]))))
        # get all values (no need to specify the id's)
        self.assertTrue(np.all(np.equal([10,20,30], self.station_manager.get_x())))
        # set a single value and get it
        self.station_manager.set_x(8, 0)
        self.assertTrue(np.all(np.equal([8,20], self.station_manager.get_x([0,1]))))
        # set two values and then get all values
        self.station_manager.set_x([16,32], [1,2])
        self.assertTrue(np.all(np.equal([8,16,32], self.station_manager.get_x())))
        
    def test_y(self):
        # get a single value from the original array
        self.assertEqual(15, self.station_manager.get_y(0))
        # get two specific values
        self.assertTrue(np.all(np.equal([25,35], self.station_manager.get_y([1,2]))))
        # get values in reverse order
        self.assertTrue(np.all(np.equal([35,25,15], self.station_manager.get_y([2,1,0]))))
        # get all values (no need to specify the id's)
        self.assertTrue(np.all(np.equal([15,25,35], self.station_manager.get_y())))
        # set a single value and get it
        self.station_manager.set_y(9, 1)
        self.assertTrue(np.all(np.equal([15,9], self.station_manager.get_y([0,1]))))
        # set two values and then get all values
        self.station_manager.set_y([7,21], [0,2])
        self.assertTrue(np.all(np.equal([7,9,21], self.station_manager.get_y())))
       
    def test_height(self):
        # get a single value from the original array
        self.assertEqual(1, self.station_manager.get_height(0))
        # get two specific values
        self.assertTrue(np.all(np.equal([1,3], self.station_manager.get_height([0,2]))))
        # get values in reverse order
        self.assertTrue(np.all(np.equal([3,2,1], self.station_manager.get_height([2,1,0]))))
        # get all values (no need to specify the id's)
        self.assertTrue(np.all(np.equal([1,2,3], self.station_manager.get_height())))
        # set a single value and get it
        self.station_manager.set_height(7, 1)
        self.assertTrue(np.all(np.equal([7,3], self.station_manager.get_height([1,2]))))
        # set two values and then get all values
        self.station_manager.set_height([5,4], [0,2])
        self.assertTrue(np.all(np.equal([5,7,4], self.station_manager.get_height())))
        
    def test_tx_power(self):
        # get a single value from the original array
        self.assertEqual(35, self.station_manager.get_tx_power(1))
        # get two specific values
        self.assertTrue(np.all(np.equal([30,40], self.station_manager.get_tx_power([0,2]))))
        # get values in reverse order
        self.assertTrue(np.all(np.equal([40,35,30], self.station_manager.get_tx_power([2,1,0]))))
        # get all values (no need to specify the id's)
        self.assertTrue(np.all(np.equal([30,35,40], self.station_manager.get_tx_power())))
        # set a single value and get it
        self.station_manager.set_tx_power(50, 1)
        self.assertTrue(np.all(np.equal([40,50], self.station_manager.get_tx_power([2,1]))))
        # set two values and then get all values
        self.station_manager.set_tx_power([20,38], [0,2])
        self.assertTrue(np.all(np.equal([20,50,38], self.station_manager.get_tx_power())))
        
    def test_rx_power(self):
        # get a single value from the original array
        self.assertEqual(-10, self.station_manager.get_rx_power(2))
        # get two specific values
        self.assertTrue(np.all(np.equal([-50,-35], self.station_manager.get_rx_power([0,1]))))
        # get values in reverse order
        self.assertTrue(np.all(np.equal([-10,-35,-50], self.station_manager.get_rx_power([2,1,0]))))
        # get all values (no need to specify the id's)
        self.assertTrue(np.all(np.equal([-50,-35,-10], self.station_manager.get_rx_power())))
        # set a single value and get it
        self.station_manager.set_rx_power(-15, 2)
        self.assertTrue(np.all(np.equal([-15,-50], self.station_manager.get_rx_power([2,0]))))
        # set two values and then get all values
        self.station_manager.set_rx_power([-60,-30], [0,1])
        self.assertTrue(np.all(np.equal([-60,-30,-15], self.station_manager.get_rx_power())))      
        
    def test_tx_antenna(self):
        self.assertEqual(10, self.station_manager.get_tx_antenna(0))
        self.assertEqual(25, self.station_manager.get_tx_antenna(1))
        self.assertEqual(Antenna(30), self.station_manager.get_tx_antenna(2))
        self.assertTrue(np.all(np.equal([10,25], self.station_manager.get_tx_antenna([0,1]))))
        self.assertTrue(np.all(np.equal([10,30], self.station_manager.get_tx_antenna([0,2]))))
        self.assertTrue(np.all(np.equal([30,25], self.station_manager.get_tx_antenna([2,1])))  )
        self.assertTrue(np.all(np.equal([30,25,10], self.station_manager.get_tx_antenna([2,1,0]))))
        self.assertTrue(np.all(np.equal([10,25,30], self.station_manager.get_tx_antenna())))
        self.assertTrue(np.all(np.equal([Antenna(10),Antenna(25),30], self.station_manager.get_tx_antenna())))
        # test setting one antenna
        self.station_manager.set_tx_antenna(Antenna(2),0)
        self.assertTrue(np.all(np.equal([2,25,30], self.station_manager.get_tx_antenna())))
        # test setting two antennas
        self.station_manager.set_tx_antenna([Antenna(2),Antenna(4)],[1,0])
        self.assertTrue(np.all(np.equal([4,2,30], self.station_manager.get_tx_antenna())))
        
    def test_rx_antenna(self):
        self.assertEqual(5, self.station_manager.get_rx_antenna(0))
        self.assertEqual(15, self.station_manager.get_rx_antenna(1))
        self.assertEqual(Antenna(20), self.station_manager.get_rx_antenna(2))
        self.assertTrue(np.all(np.equal([5,15], self.station_manager.get_rx_antenna([0,1]))))
        self.assertTrue(np.all(np.equal([5,20], self.station_manager.get_rx_antenna([0,2]))))
        self.assertTrue(np.all(np.equal([20,15], self.station_manager.get_rx_antenna([2,1]))))
        self.assertTrue(np.all(np.equal([20,15,5], self.station_manager.get_rx_antenna([2,1,0]))))
        self.assertTrue(np.all(np.equal([5,15,20], self.station_manager.get_rx_antenna())))
        self.assertTrue(np.all(np.equal([Antenna(5),Antenna(15),20], self.station_manager.get_rx_antenna())))
        # test setting one antenna
        self.station_manager.set_rx_antenna(Antenna(7),0)
        self.assertTrue(np.all(np.equal([7,15,20], self.station_manager.get_rx_antenna())))
        # test setting two antennas
        self.station_manager.set_rx_antenna([Antenna(8),Antenna(9)],[1,0])
        self.assertTrue(np.all(np.equal([9,8,20], self.station_manager.get_rx_antenna())))
        
    def test_station(self):
        # test if manager returns the correct station
        self.assertTrue(self.station == self.station_manager.get_station(0))
        self.assertTrue(self.station2 == self.station_manager.get_station(1))
        # now we change station id and verify if stations became different
        self.station.set_id(1)
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
        ref_distance = np.asarray([[ 356.405,  180.277]])
        distance = self.station_manager3.get_distance_to(self.station_manager2)
        self.assertTrue(np.all(np.isclose(ref_distance, distance, atol=1e-2)))        
        
        ref_distance = np.asarray([[ 127.279,  302.200],
                                   [ 113.137,  288.140],
                                   [  98.994,  274.089]])
        distance = self.station_manager.get_distance_to(self.station_manager2)
        self.assertTrue(np.all(np.isclose(ref_distance, distance, atol=1e-2)))
        
    def test_3d_distance_to(self):
        ref_distance = np.asarray([[ 356.411,  180.302]])
        distance = self.station_manager3.get_3d_distance_to(self.station_manager2)
        self.assertTrue(np.all(np.isclose(ref_distance, distance, atol=1e-2)))        
        
        ref_distance = np.asarray([[ 127.314,  302.226],
                                   [ 113.154,  288.156],
                                   [  99,  274.096]])
        distance = self.station_manager.get_3d_distance_to(self.station_manager2)
        self.assertTrue(np.all(np.isclose(ref_distance, distance, atol=1e-2)))        
        
if __name__ == '__main__':
    unittest.main()
                