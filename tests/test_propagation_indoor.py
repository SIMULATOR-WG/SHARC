# -*- coding: utf-8 -*-
"""
Created on Mon Dec 18 20:10:14 2017

@author: edgar
"""

import unittest
import numpy as np
import numpy.testing as npt

from sharc.propagation.propagation_indoor import PropagationIndoor
from sharc.parameters.parameters_indoor import ParametersIndoor


class PropagationIndoorTest(unittest.TestCase):

    def setUp(self):
        pass

    def test_loss(self):
        params = ParametersIndoor()
        params.basic_path_loss = "INH_OFFICE"
        params.n_rows = 3
        params.n_colums = 1
    #    params.street_width = 30
        params.ue_indoor_percent = .95
        params.building_class = "TRADITIONAL"
        params.num_cells = 3

        bs_per_building = 3
        ue_per_bs = 3

        num_bs = bs_per_building*params.n_rows*params.n_colums
        num_ue = num_bs*ue_per_bs
        distance_2D = 150*np.random.random((num_bs, num_ue))
        frequency = 27000*np.ones(distance_2D.shape)
        indoor = np.random.rand(num_bs) < params.ue_indoor_percent
        h_bs = 3*np.ones(num_bs)
        h_ue = 1.5*np.ones(num_ue)
        distance_3D = np.sqrt(distance_2D**2 + (h_bs[:,np.newaxis] - h_ue)**2)
        height_diff = np.tile(h_bs, (num_bs, 3)) - np.tile(h_ue, (num_bs, 1))
        elevation = np.degrees(np.arctan(height_diff/distance_2D))

        propagation_indoor = PropagationIndoor(np.random.RandomState(), params, ue_per_bs)

if __name__ == '__main__':
    unittest.main()
