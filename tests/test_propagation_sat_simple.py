# -*- coding: utf-8 -*-
"""
Created on Wed Jul 12 13:56:47 2017

@author: edgar
"""

import unittest
import numpy.testing as npt

from sharc.propagation.propagation_sat_simple import PropagationSatSimple

class PropagationSatSimpleTest(unittest.TestCase):
    
    
    def setUp(self):
        self.propagation = PropagationSatSimple()
        
        
    def test_loss_los(self):
        d = 10
        f = 10
        los_prob = 1
        self.assertEqual(12.45 + 4, 
                         self.propagation.get_loss(distance_3D=d, 
                                                   frequency=f,
                                                   line_of_sight_prob=los_prob))

        d = [ 10, 100 ]
        f = [ 10, 100 ]
        npt.assert_allclose([12.45 + 4, 52.45 + 4], 
                            self.propagation.get_loss(distance_3D=d, 
                                                      frequency=f,
                                                      line_of_sight_prob=los_prob),
                            atol=1e-2)

        d = [ 10, 100, 1000 ]
        f = [ 10, 100, 1000 ]
        npt.assert_allclose([12.45 + 4, 52.45 + 4, 92.45 + 4], 
                            self.propagation.get_loss(distance_3D=d, 
                                                      frequency=f,
                                                      line_of_sight_prob=los_prob),
                            atol=1e-2)

        d = [[10, 20, 30],[40, 50, 60]]
        f = [ 100 ]
        ref_loss = [[ 32.45 + 4,  38.47 + 4,  41.99 + 4],
                    [ 44.49 + 4,  46.42 + 4,  48.01 + 4]]
        loss = self.propagation.get_loss(distance_3D=d, 
                                         frequency=f,
                                         line_of_sight_prob=los_prob)
        npt.assert_allclose(ref_loss, loss, atol=1e-2)   
        
        
    def test_loss_nlos(self):
        d = 10
        f = 10
        los_prob = 0
        self.assertEqual(12.45 + 4 + 20, 
                         self.propagation.get_loss(distance_3D=d, 
                                                   frequency=f,
                                                   line_of_sight_prob=los_prob))

        d = [ 10, 100 ]
        f = [ 10, 100 ]
        npt.assert_allclose([12.45 + 4 + 20, 52.45 + 4 + 20], 
                            self.propagation.get_loss(distance_3D=d, 
                                                      frequency=f,
                                                      line_of_sight_prob=los_prob),
                            atol=1e-2)

        d = [ 10, 100, 1000 ]
        f = [ 10, 100, 1000 ]
        npt.assert_allclose([12.45 + 4 + 20, 52.45 + 4 + 20, 92.45 + 4 + 20], 
                            self.propagation.get_loss(distance_3D=d, 
                                                      frequency=f,
                                                      line_of_sight_prob=los_prob),
                            atol=1e-2)

        d = [[10, 20, 30],[40, 50, 60]]
        f = [ 100 ]
        ref_loss = [[ 32.45 + 4 + 20,  38.47 + 4 + 20,  41.99 + 4 + 20],
                    [ 44.49 + 4 + 20,  46.42 + 4 + 20,  48.01 + 4 + 20]]
        loss = self.propagation.get_loss(distance_3D=d, 
                                         frequency=f,
                                         line_of_sight_prob=los_prob)
        npt.assert_allclose(ref_loss, loss, atol=1e-2)         
        
if __name__ == '__main__':
    unittest.main()
        