# -*- coding: utf-8 -*-
"""
Created on Mon Feb  6 17:59:07 2017

@author: edgar
"""

import unittest

from antenna import Antenna
from station import Station

class StationTest(unittest.TestCase):
    
    def setUp(self):
        self.station = Station()
        self.station.set_id(1)
        self.station.set_x(10)
        self.station.set_y(15)
        self.station.set_height(6)
        self.station.set_tx_power(20)
        self.station.set_rx_power(-3)
        self.station.set_tx_antenna(Antenna(30))
        self.station.set_rx_antenna(Antenna(35))
        
        self.station2 = Station()
        self.station2.set_id(1)
        self.station2.set_x(10)
        self.station2.set_y(15)
        self.station2.set_height(6)
        self.station2.set_tx_power(17)
        self.station2.set_rx_power(9)
        self.station2.set_tx_antenna(Antenna(10))
        self.station2.set_rx_antenna(Antenna(12))       
        
        self.station3 = Station()
        self.station3.set_id(2)
        self.station3.set_x(10)
        self.station3.set_y(15)
        self.station3.set_height(6)
        self.station3.set_tx_power(20)
        self.station3.set_rx_power(-3)
        self.station3.set_tx_antenna(Antenna(30))
        self.station3.set_rx_antenna(Antenna(35))       
        
    def test_id(self):
        self.assertEqual(1, self.station.get_id())

    def test_x(self):
        self.assertEqual(10, self.station.get_x())

    def test_y(self):
        self.assertEqual(15, self.station.get_y())

    def test_height(self):
        self.assertEqual(6, self.station.get_height())

    def test_tx_power(self):
        self.assertEqual(20, self.station.get_tx_power())

    def test_rx_power(self):
        self.assertEqual(-3, self.station.get_rx_power())

    def test_tx_antenna(self):
        self.assertEqual(30, self.station.get_tx_antenna().get_gain())

    def test_rx_antenna(self):
        self.assertEqual(35, self.station.get_rx_antenna().get_gain())
        
    def test_eq(self):
        self.assertTrue(self.station == self.station2)
        # changing id, x, y, or height should change the result
        self.station.set_x(11)
        self.assertFalse(self.station == self.station2)
        #
        self.assertFalse(self.station == self.station3)
        self.assertFalse(self.station2 == self.station3)
        
    def test_ne(self):
        self.assertFalse(self.station != self.station2)
        # changing id, x, y, or height should change the result
        self.station.set_height(9)
        self.assertTrue(self.station != self.station2)   
        #
        self.assertTrue(self.station != self.station3)
        self.assertTrue(self.station2 != self.station3)
        
        
        
if __name__ == '__main__':
    unittest.main()
        