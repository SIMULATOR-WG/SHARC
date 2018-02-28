# -*- coding: utf-8 -*-
"""
Created on Mon Aug 21 12:17:34 2017

@author: Andre Barreto
"""

import unittest
import numpy as np
import numpy.testing as npt
from sharc.propagation.propagation_building_entry_loss import PropagationBuildingEntryLoss


class TestPropagationBuildingEntryLoss(unittest.TestCase):

    def setUp(self):
        self.building_entry_loss = PropagationBuildingEntryLoss(np.random.RandomState())

    def test_building_entry_loss(self):
        # compare with benchmark from ITU-R P-2109-0 Fig. 1

        f_GHz_vec = np.array([.1, .3, 1, 5, 10, 30, 100])

        # traditional building type
        loss_lower = np.array([9, 10, 11, 15, 16, 19, 21])
        loss_upper = np.array([10, 12, 15, 17, 18, 21, 25])


        loss = self.building_entry_loss.get_loss(f_GHz_vec * 1000, 0, prob=.5, test=True)
        npt.assert_array_less(loss_lower, loss)
        npt.assert_array_less(loss, loss_upper)

        # energy-efficient building type
        loss_lower = [39, 30, 25, 30, 32, 40, 55]
        loss_upper = [40, 35, 30, 31, 35, 43, 57]

        loss = self.building_entry_loss.get_loss(f_GHz_vec * 1000, 0, prob=.5, test=True,
                                                 building_class="THERMALLY_EFFICIENT")
        npt.assert_array_less(loss_lower, loss)
        npt.assert_array_less(loss, loss_upper)


if __name__ == '__main__':
    unittest.main()
