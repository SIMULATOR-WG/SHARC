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

        d = [[10, 20, 30],[40, 50, 60]]
        f = [ 100 ]
        ref_loss = [[ 32.45,  38.47,  41.99],
                    [ 44.49,  46.42,  48.01]]
        loss = self.__freeSpace.get_loss(distance=d, frequency=f)
        self.assertTrue(np.all(np.isclose(ref_loss, loss, atol=1e-2)))         
        
        
if __name__ == '__main__':
    unittest.main()
        