# -*- coding: utf-8 -*-
"""
Created on Fri Aug 10 17:38:33 2018

@author: Calil
"""

import unittest
import numpy as np
import numpy.testing as npt

from sharc.parameters.parameters_fss_es import ParametersFssEs
from sharc.support.enumerations import StationType
from sharc.propagation.propagation_hdfss import PropagationHDFSS

class PropagationHDFSSTest(unittest.TestCase):

    def setUp(self):
        rnd = np.random.RandomState(101)
        par = ParametersFssEs()
        par.building_loss_enabled = False
        par.bs_building_entry_loss_type = 'FIXED_VALUE'
        par.bs_building_entry_loss_prob = 0.5
        par.bs_building_entry_loss_value = 50
        self.propagation = PropagationHDFSS(par,rnd)
        
        # Propagation with fixed BEL
        rnd = np.random.RandomState(101)
        par = ParametersFssEs()
        par.building_loss_enabled = True
        par.bs_building_entry_loss_type = 'FIXED_VALUE'
        par.bs_building_entry_loss_prob = 0.6
        par.bs_building_entry_loss_value = 50
        self.propagation_fixed_value = PropagationHDFSS(par,rnd)
        
        # Propagation with fixed probability
        rnd = np.random.RandomState(101)
        par = ParametersFssEs()
        par.building_loss_enabled = True
        par.bs_building_entry_loss_type = 'P2109_FIXED'
        par.bs_building_entry_loss_prob = 0.6
        par.bs_building_entry_loss_value = 50
        self.propagation_fixed_prob = PropagationHDFSS(par,rnd)
        
        # Propagation with random probability
        rnd = np.random.RandomState(101)
        par = ParametersFssEs()
        par.building_loss_enabled = True
        par.bs_building_entry_loss_type = 'P2109_RANDOM'
        par.bs_building_entry_loss_prob = 0.6
        par.bs_building_entry_loss_value = 50
        self.propagation_random_prob = PropagationHDFSS(par,rnd)
        
    def test_get_loss(self):
        d = np.array([10.0, 20.0, 30.0, 60.0, 90.0, 300.0, 1000.0])
        f = 40000*np.ones_like(d)
        ele = np.zeros_like(d)
        
        loss = self.propagation.get_loss(distance_3D=d,
                                         frequency=f,
                                         elevation=ele,
                                         imt_sta_type=StationType.IMT_BS,
                                         shadow=False)
        
        expected_loss = np.array([84.48, 90.50, 94.02, 100.72, 104.75, 139.33, 162.28])
        
        npt.assert_allclose(loss,expected_loss,atol=1e-1)
    
    def test_get_build_loss(self):
        # Initialize variables
        ele = np.array([ 0.0, 45.0, 90.0])
        f = 40000*np.ones_like(ele)
        sta_type = StationType.IMT_BS
        
        # Test 1: fixed value
        expected_build_loss = 50.0
        build_loss = self.propagation_fixed_value.get_building_loss(sta_type,
                                                                    f,
                                                                    ele)
        self.assertEqual(build_loss,expected_build_loss)
        
        # Test 2: fixed probability
        expected_build_loss = np.array([24.4, 33.9, 43.4])
        build_loss = self.propagation_fixed_prob.get_building_loss(sta_type,
                                                                   f,
                                                                   ele)
        npt.assert_allclose(build_loss,expected_build_loss,atol=1e-1)
        
        # Test 3: random probability
        expected_build_loss = np.array([21.7, 32.9, 15.9])
        build_loss = self.propagation_random_prob.get_building_loss(sta_type,
                                                                    f,
                                                                    ele)
        npt.assert_allclose(build_loss,expected_build_loss,atol=1e-1)
        
        # Test 4: UE station
        sta_type = StationType.IMT_UE
        expected_build_loss = np.array([21.7, 32.9, 15.9])
        build_loss = self.propagation_fixed_value.get_building_loss(sta_type,
                                                                    f,
                                                                    ele)
        npt.assert_allclose(build_loss,expected_build_loss,atol=1e-1)
        build_loss = self.propagation_fixed_prob.get_building_loss(sta_type,
                                                                    f,
                                                                    ele)
        npt.assert_allclose(build_loss,expected_build_loss,atol=1e-1)
        expected_build_loss = np.array([10.1, 36.8, 52.6])
        build_loss = self.propagation_random_prob.get_building_loss(sta_type,
                                                                    f,
                                                                    ele)
        npt.assert_allclose(build_loss,expected_build_loss,atol=1e-1)
        
    
if __name__ == '__main__':
    unittest.main()
