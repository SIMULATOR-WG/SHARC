# -*- coding: utf-8 -*-
"""
Created on Wed Jul 12 13:56:47 2017

@author: edgar
"""

import unittest
import numpy.testing as npt
import numpy as np

from sharc.propagation.propagation_sat_simple import PropagationSatSimple

class PropagationSatSimpleTest(unittest.TestCase):
    
    
    def setUp(self):
        self.propagation = PropagationSatSimple()
        
        
    def test_loss_los(self):
        d = np.array(10)
        f = np.array(10)
        los_prob = 1
        indoor_stations = np.zeros(d.shape, dtype=bool)
        self.assertEqual(12.45 + 4, 
                         self.propagation.get_loss(distance_3D=d, 
                                                   frequency=f,
                                                   indoor_stations=indoor_stations,
                                                   line_of_sight_prob=los_prob))

        d = np.array([ 10, 100 ])
        f = np.array([ 10, 100 ])
        indoor_stations = np.zeros(d.shape, dtype=bool)
        npt.assert_allclose([12.45 + 4, 52.45 + 4], 
                            self.propagation.get_loss(distance_3D=d, 
                                                      frequency=f,
                                                      indoor_stations=indoor_stations,
                                                      line_of_sight_prob=los_prob),
                            atol=1e-2)

        d = np.array([ 10, 100, 1000 ])
        f = np.array([ 10, 100, 1000 ])
        indoor_stations = np.array([1, 0, 0], dtype=bool)
        npt.assert_allclose([12.45 + 4 + 20, 52.45 + 4, 92.45 + 4], 
                            self.propagation.get_loss(distance_3D=d, 
                                                      frequency=f,
                                                      indoor_stations=indoor_stations,
                                                      line_of_sight_prob=los_prob),
                            atol=1e-2)

        d = np.array([[10, 20, 30],[40, 50, 60]])
        f = np.array([ 100 ])
        indoor_stations = np.zeros(d.shape, dtype=bool)
        ref_loss = [[ 32.45 + 4,  38.47 + 4,  41.99 + 4],
                    [ 44.49 + 4,  46.42 + 4,  48.01 + 4]]
        loss = self.propagation.get_loss(distance_3D=d, 
                                         frequency=f,
                                         indoor_stations=indoor_stations,
                                         line_of_sight_prob=los_prob)
        npt.assert_allclose(ref_loss, loss, atol=1e-2)   
        
        
    def test_loss_nlos(self):
        d = np.array(10)
        f = np.array(10)
        los_prob = 0
        indoor_stations = np.zeros(d.shape, dtype=bool)
        self.assertEqual(12.45 + 4 + 20, 
                         self.propagation.get_loss(distance_3D=d, 
                                                   frequency=f,
                                                   indoor_stations=indoor_stations,
                                                   line_of_sight_prob=los_prob))

        d = np.array([ 10, 100 ])
        f = np.array([ 10, 100 ])
        indoor_stations = np.zeros(d.shape, dtype=bool)
        npt.assert_allclose([12.45 + 4 + 20, 52.45 + 4 + 20], 
                            self.propagation.get_loss(distance_3D=d, 
                                                      frequency=f,
                                                      indoor_stations=indoor_stations,
                                                      line_of_sight_prob=los_prob),
                            atol=1e-2)

        d = np.array([ 10, 100, 1000 ])
        f = np.array([ 10, 100, 1000 ])
        indoor_stations = np.array([0, 0, 1], dtype=bool)
        npt.assert_allclose([12.45 + 4 + 20, 52.45 + 4 + 20, 92.45 + 4 + 20 + 20], 
                            self.propagation.get_loss(distance_3D=d, 
                                                      frequency=f,
                                                      indoor_stations=indoor_stations,
                                                      line_of_sight_prob=los_prob),
                            atol=1e-2)

        d = np.array([[10, 20, 30],[40, 50, 60]])
        f = np.array([ 100 ])
        indoor_stations = np.zeros(d.shape, dtype=bool)
        ref_loss = [[ 32.45 + 4 + 20,  38.47 + 4 + 20,  41.99 + 4 + 20],
                    [ 44.49 + 4 + 20,  46.42 + 4 + 20,  48.01 + 4 + 20]]
        loss = self.propagation.get_loss(distance_3D=d, 
                                         frequency=f,
                                         indoor_stations=indoor_stations,
                                         line_of_sight_prob=los_prob)
        npt.assert_allclose(ref_loss, loss, atol=1e-2)         
        
if __name__ == '__main__':
    unittest.main()
        