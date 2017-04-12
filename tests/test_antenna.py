# -*- coding: utf-8 -*-
"""
Created on Mon Feb  6 16:39:34 2017

@author: edgar
"""

import unittest

from sharc.antenna import Antenna

class AntennaTest(unittest.TestCase):
    
    def setUp(self):
        self.antenna1 = Antenna()
        self.antenna1.gain = 5
        
        self.antenna2 = Antenna()
        self.antenna2.gain = 8.0
        
        self.antenna3 = Antenna(10)
        
    def test_gain(self):
        self.assertEqual(self.antenna1.gain, 5)
        self.assertEqual(self.antenna2.gain, 8)
        self.assertEqual(self.antenna3.gain, 10)
        
    def test_float(self):
        self.assertTrue(type(float(self.antenna1)) is float )
        self.assertTrue(type(float(self.antenna2)) is float )
        
    def test_add(self):
        self.assertEqual(self.antenna1 + 2, 7)
        self.assertEqual(self.antenna2 + 3, 11)
    
    def test_radd(self):
        self.assertEqual(2 + self.antenna1, 7)
        self.assertEqual(3 + self.antenna2, 11)
        self.assertEqual(self.antenna1 + self.antenna2, 13)
        
    def test_sub(self):
        self.assertEqual(self.antenna1 - 1, 4)
        self.assertEqual(self.antenna2 - 2, 6)
    
    def test_rsub(self):
        self.assertEqual(9 - self.antenna1, 4)
        self.assertEqual(10 - self.antenna2, 2)
     
    def test_lt(self):
        self.assertFalse(self.antenna1 < 5)
        self.assertTrue(self.antenna2 < 10)
        self.assertTrue(self.antenna2 < self.antenna3)
        
    def test_le(self):
        self.assertTrue(self.antenna1 <= 5)
        self.assertFalse(self.antenna2 <= 7.9)
        self.assertTrue(self.antenna1 <= self.antenna3)

    def test_gt(self):
        self.assertFalse(self.antenna1 > 5)
        self.assertTrue(self.antenna2 > 7)
        self.assertTrue(self.antenna3 > self.antenna1)
        
    def test_ge(self):
        self.assertTrue(self.antenna1 >= 5)
        self.assertFalse(self.antenna2 >= 8.1)
        self.assertTrue(self.antenna3 >= self.antenna1)
        
    def test_eq(self):
        self.assertTrue(self.antenna1 == 5)
        self.assertFalse(self.antenna2 == 7)
        self.assertTrue(self.antenna3 == Antenna(10))
        
    def test_ne(self):
        self.assertTrue(self.antenna1 != 1)
        self.assertFalse(self.antenna2 != 8)
        self.assertTrue(self.antenna3 != Antenna(3))        
        
        
if __name__ == '__main__':
    unittest.main()
