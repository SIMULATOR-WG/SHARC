# -*- coding: utf-8 -*-
"""
Created on Mon Mar  13 15:14:34 2017

@author: edgar
"""

import unittest
import numpy as np

from propagation_free_space import PropagationFreeSpace

class PropagationFreeSpaceTest(unittest.TestCase):
    
    def setUp(self):
        self.__freeSpace = PropagationFreeSpace()
        
    def test_loss(self):
        d = 10
        f = 10
        self.assertEqual(12.45, 
                         self.__freeSpace.get_loss(distance=d, frequency=f))

        d = [ 10, 100 ]
        f = [ 10, 100 ]
        self.assertTrue(np.all(np.equal([12.45, 52.45], 
                         self.__freeSpace.get_loss(distance=d, frequency=f))))

        d = [ 10, 100, 1000 ]
        f = [ 10, 100, 1000 ]
        self.assertTrue(np.all(np.equal([12.45, 52.45, 92.45], 
                         self.__freeSpace.get_loss(distance=d, frequency=f))))        
        
        
if __name__ == '__main__':
    unittest.main()
        