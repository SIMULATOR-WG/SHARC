# -*- coding: utf-8 -*-
"""
Created on Wed Aug 29 16:27:33 2018

@author: Calil
"""

import unittest
import numpy as np
import numpy.testing as npt

from sharc.parameters.parameters_fss_es import ParametersFssEs
from sharc.support.enumerations import StationType
from sharc.propagation.propagation_hdfss_building_side import PropagationHDFSSBuildingSide

class PropagationHDFSSBuildingSideTest(unittest.TestCase):

    def setUp(self):
        # Basic propagation
        rnd = np.random.RandomState(101)
        par = ParametersFssEs()
        par.building_loss_enabled = False
        par.shadow_enabled = False
        par.same_building_enabled = True
        par.diffraction_enabled = False
        par.bs_building_entry_loss_type = 'FIXED_VALUE'
        par.bs_building_entry_loss_prob = 0.5
        par.bs_building_entry_loss_value = 50
        self.propagation = PropagationHDFSSBuildingSide(par,rnd)
        
        # Propagation with fixed BEL
        rnd = np.random.RandomState(101)
        par = ParametersFssEs()
        par.building_loss_enabled = True
        par.shadow_enabled = False
        par.same_building_enabled = True
        par.diffraction_enabled = False
        par.bs_building_entry_loss_type = 'FIXED_VALUE'
        par.bs_building_entry_loss_prob = 0.6
        par.bs_building_entry_loss_value = 50
        self.propagation_fixed_value = PropagationHDFSSBuildingSide(par,rnd)
        
        # Propagation with fixed probability
        rnd = np.random.RandomState(101)
        par = ParametersFssEs()
        par.building_loss_enabled = True
        par.shadow_enabled = False
        par.same_building_enabled = True
        par.diffraction_enabled = False
        par.bs_building_entry_loss_type = 'P2109_FIXED'
        par.bs_building_entry_loss_prob = 0.6
        par.bs_building_entry_loss_value = 50
        self.propagation_fixed_prob = PropagationHDFSSBuildingSide(par,rnd)
        
        # Propagation with random probability
        rnd = np.random.RandomState(101)
        par = ParametersFssEs()
        par.building_loss_enabled = True
        par.shadow_enabled = False
        par.same_building_enabled = True
        par.diffraction_enabled = False
        par.bs_building_entry_loss_type = 'P2109_RANDOM'
        par.bs_building_entry_loss_prob = 0.6
        par.bs_building_entry_loss_value = 50
        self.propagation_random_prob = PropagationHDFSSBuildingSide(par,rnd)
        
    def test_get_loss(self):
        pass
    
    def test_get_build_loss(self):
        # Initialize variables
        ele = np.array([[ 0.0, 45.0, 90.0]])
        f = 40000*np.ones_like(ele)
        sta_type = StationType.IMT_BS
        
        # Test 1: fixed value
        expected_build_loss = 50.0
        build_loss = self.propagation_fixed_value.get_building_loss(sta_type,
                                                                    f,
                                                                    ele)
        self.assertEqual(build_loss,expected_build_loss)
        
        # Test 2: fixed probability
        expected_build_loss = np.array([[24.4, 33.9, 43.4]])
        build_loss = self.propagation_fixed_prob.get_building_loss(sta_type,
                                                                   f,
                                                                   ele)
        npt.assert_allclose(build_loss,expected_build_loss,atol=1e-1)
        
        # Test 3: random probability
        expected_build_loss = np.array([[21.7, 32.9, 15.9]])
        build_loss = self.propagation_random_prob.get_building_loss(sta_type,
                                                                    f,
                                                                    ele)
        npt.assert_allclose(build_loss,expected_build_loss,atol=1e-1)
        
        # Test 4: UE station
        sta_type = StationType.IMT_UE
        expected_build_loss = np.array([[21.7, 32.9, 15.9]])
        build_loss = self.propagation_fixed_value.get_building_loss(sta_type,
                                                                    f,
                                                                    ele)
        npt.assert_allclose(build_loss,expected_build_loss,atol=1e-1)
        build_loss = self.propagation_fixed_prob.get_building_loss(sta_type,
                                                                    f,
                                                                    ele)
        npt.assert_allclose(build_loss,expected_build_loss,atol=1e-1)
        expected_build_loss = np.array([[10.1, 36.8, 52.6]])
        build_loss = self.propagation_random_prob.get_building_loss(sta_type,
                                                                    f,
                                                                    ele)
        npt.assert_allclose(build_loss,expected_build_loss,atol=1e-1)
    
    def test_is_same_build(self):
        # Test is_same_building()
        es_x = np.array([0.0])
        es_y = np.array([25.0])
        imt_x = np.array([1.0, 0.0,80.0,-70.0,12.0])
        imt_y = np.array([1.0,30.0, 0.0,-29.3,-3.6])
        
        expected_in_build = np.array([True,False,False,False,True])
        in_build = self.propagation.is_same_building(imt_x,
                                                     imt_y,
                                                     es_x,
                                                     es_y)
        npt.assert_array_equal(in_build,expected_in_build)
    
    def test_is_next_build(self):
        # Test is_same_building()
        es_x = np.array([0.0])
        es_y = np.array([25.0])
        imt_x = np.array([1.0, 0.0,80.0,-70.0,12.0])
        imt_y = np.array([1.0,30.0, 0.0,-29.3,-3.6]) + 80.0
        
        expected_in_build = np.array([True,False,False,False,True])
        in_build = self.propagation.is_next_building(imt_x,
                                                     imt_y,
                                                     es_x,
                                                     es_y)
        npt.assert_array_equal(in_build,expected_in_build)
        
    
if __name__ == '__main__':
    unittest.main()