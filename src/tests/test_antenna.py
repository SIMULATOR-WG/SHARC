# -*- coding: utf-8 -*-
"""
Created on Mon Feb  6 16:39:34 2017

@author: edgar
"""

import unittest

from antenna import Antenna

class AntennaTest(unittest.TestCase):
    
    def setUp(self):
        self.antenna1 = Antenna()
        self.antenna1.set_gain(5)
        self.antenna2 = Antenna()
        self.antenna2.set_gain(8.0)
        self.antenna3 = Antenna(10)
        
    def test_gain(self):
        self.assertEqual(5, self.antenna1.get_gain())
        self.assertEqual(8, self.antenna2.get_gain())
        self.assertEqual(10, self.antenna3.get_gain())
        
    def test_float(self):
        self.assertTrue(type(float(self.antenna1)) is float )
        self.assertTrue(type(float(self.antenna2)) is float )
        
    def test_add(self):
        self.assertEqual(7, self.antenna1 + 2)
        self.assertEqual(11, self.antenna2 + 3)
    
    def test_radd(self):
        self.assertEqual(7, 2 + self.antenna1)
        self.assertEqual(11, 3 + self.antenna2)
        self.assertEqual(13, self.antenna1 + self.antenna2)
        
    def test_sub(self):
        self.assertEqual(4, self.antenna1 - 1)
        self.assertEqual(6, self.antenna2 - 2)
    
    def test_rsub(self):
        self.assertEqual(4, 9 - self.antenna1)
        self.assertEqual(2, 10 - self.antenna2)
     
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
