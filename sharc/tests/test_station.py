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
        self.station.id = 1
        self.station.x = 10
        self.station.y = 15
        self.station.height = 6
        self.station.tx_power = 20
        self.station.rx_power = -3
        self.station.tx_antenna = Antenna(30)
        self.station.rx_antenna = Antenna(35)
        
        self.station2 = Station()
        self.station2.id = 1
        self.station2.x = 10
        self.station2.y = 15
        self.station2.height = 6
        self.station2.tx_power = 17
        self.station2.rx_power = 9
        self.station2.tx_antenna = Antenna(10)
        self.station2.rx_antenna = Antenna(12)
        
        self.station3 = Station()
        self.station3.id = 2
        self.station3.x = 10
        self.station3.y = 15
        self.station3.height = 6
        self.station3.tx_power = 20
        self.station3.rx_power = -3
        self.station3.tx_antenna = Antenna(30)
        self.station3.rx_antenna = Antenna(35)   
        
    def test_id(self):
        self.assertEqual(self.station.id, 1)

    def test_x(self):
        self.assertEqual(self.station.x, 10)

    def test_y(self):
        self.assertEqual(self.station.y, 15)

    def test_height(self):
        self.assertEqual(self.station.height, 6)

    def test_tx_power(self):
        self.assertEqual(self.station.tx_power, 20)

    def test_rx_power(self):
        self.assertEqual(self.station.rx_power, -3)

    def test_tx_antenna(self):
        self.assertEqual(self.station.tx_antenna.gain, 30)

    def test_rx_antenna(self):
        self.assertEqual(self.station.rx_antenna.gain, 35)
        
    def test_eq(self):
        self.assertTrue(self.station == self.station2)
        # changing id, x, y, or height should change the result
        self.station.x = 11
        self.assertFalse(self.station == self.station2)
        #
        self.assertFalse(self.station == self.station3)
        self.assertFalse(self.station2 == self.station3)
        
    def test_ne(self):
        self.assertFalse(self.station != self.station2)
        # changing id, x, y, or height should change the result
        self.station.height = 9
        self.assertTrue(self.station != self.station2)   
        #
        self.assertTrue(self.station != self.station3)
        self.assertTrue(self.station2 != self.station3)
        
        
        
if __name__ == '__main__':
    unittest.main()
        