# -*- coding: utf-8 -*-

"""
Created on  Mar  04 14:13:31 2017

@author:  LeticiaValle_Mac
"""

import unittest
import numpy as np
import numpy.testing as npt


from sharc.propagation.propagation_abg import PropagationABG


class PropagationABGTest(unittest.TestCase):

    def setUp(self):
        self.abg = PropagationABG(
            random_number_gen=np.random.RandomState(),
            alpha=3.4,
            beta=19.2,
            gamma=2.3,
            building_loss=20,
            shadowing_sigma_dB=6.5,
        )

    def test_loss(self):
        d = np.array([[100, 500], [400, 60]])
        f = np.ones(d.shape, dtype=float) * 27000.0
        indoor = np.zeros(d.shape[0], dtype=bool)
        shadowing = False
        loss = np.array([[120.121, 143.886347], [140.591406, 112.578509]])

        npt.assert_allclose(
            self.abg.get_loss(d, f, indoor, shadowing),
            loss, atol=1e-2,
        )

        d = np.array([500, 3000])[:, np.newaxis]
        f = np.array([27000, 40000])[:, np.newaxis]
        indoor = np.zeros(d.shape[0], dtype=bool)
        shadowing = False
        loss = np.array([143.886, 174.269])[:, np.newaxis]
        npt.assert_allclose(
            self.abg.get_loss(d, f, indoor, shadowing),
            loss, atol=1e-2,
        )


if __name__ == '__main__':
    unittest.main()
