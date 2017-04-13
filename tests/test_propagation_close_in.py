# -*- coding: utf-8 -*-
"""
Created on Mon Mar  13 15:14:34 2017

@author: edgar
"""

import unittest
import numpy as np

from sharc.propagation.propagation_close_in import PropagationCloseIn

class PropagationCloseInTest(unittest.TestCase):

    def setUp(self):
        self.__closeInMacroLos = PropagationCloseIn("MACROCELL", True, False)
        self.__closeInMacroNLos = PropagationCloseIn("MACROCELL", False, False)
        self.__closeInMacroShadowing = PropagationCloseIn("MACROCELL", True, True)

    def test_loss_macro_los(self):
        d = 10
        f = 10
        self.assertEqual(-27.55+20+20,
                         self.__closeInMacroLos.get_loss(distance=d, frequency=f))

        d = [ 10, 100 ]
        f = [ 10, 100 ]
        self.assertTrue(np.all(np.equal([-27.55+20+20, -27.55+40+40],
                         self.__closeInMacroLos.get_loss(distance=d, frequency=f))))

        d = [ 10, 100, 1000 ]
        f = [ 10, 100, 1000 ]
        self.assertTrue(np.all(np.equal([-27.55+20+20, -27.55+40+40, -27.55+60+60 ],
                         self.__closeInMacroLos.get_loss(distance=d, frequency=f))))

        d = [[1, 10],[100, 1000]]
        f = [ 100 ]
        ref_loss = [[ -27.55+0+40,  -27.55+20+40 ],
                    [ -27.55+40+40, -27.55+60+40]]
        loss = self.__closeInMacroLos.get_loss(distance=d, frequency=f)
        self.assertTrue(np.all(np.equal(ref_loss, loss)))


    def test_loss_nlos(self):
        d = 10
        f = 10
        self.assertEqual(-27.55 + 30 + 20,
                         self.__closeInMacroNLos.get_loss(distance=d, frequency=f))

        d = [10, 100]
        f = [10, 100]
        self.assertTrue(np.all(np.equal([-27.55 + 30 + 20, -27.55 + 60 + 40],
                                        self.__closeInMacroNLos.get_loss(distance=d, frequency=f))))

        d = [10, 100, 1000]
        f = [10, 100, 1000]
        self.assertTrue(np.all(np.equal([-27.55 + 30 + 20, -27.55 + 60 + 40, -27.55 + 90 + 60],
                                        self.__closeInMacroNLos.get_loss(distance=d, frequency=f))))

        d = [[1, 10], [100, 1000]]
        f = [100]
        ref_loss = [[-27.55 + 0 + 40, -27.55 + 30 + 40],
                    [-27.55 + 60 + 40, -27.55 + 90 + 40]]
        loss = self.__closeInMacroNLos.get_loss(distance=d, frequency=f)
        self.assertTrue(np.all(np.equal(ref_loss, loss)))


    def test_shadowing (self):
        d = np.ones(100000)
        f = np.ones(100000)
        loss = self.__closeInMacroShadowing.get_loss(distance=d, frequency=f)

        mean = np.mean(loss)
        std = np.std(loss)

        self.assertAlmostEqual(mean, -27.55, delta = .1)
        self.assertAlmostEqual(std, 4.1, delta = .1)


if __name__ == '__main__':
    unittest.main()
