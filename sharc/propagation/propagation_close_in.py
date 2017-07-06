# -*- coding: utf-8 -*-
"""
Created on Wed Apr 12 16:59:27 2017

@author: andre barreto
"""

from sharc.propagation.propagation import Propagation

import numpy as np

class PropagationCloseIn(Propagation):
    """
    Implements the close-in free-space reference distance propagation model.
    """

    def __init__(self):
        super().__init__()
        # model parameters for LOS
        self.path_loss_exponent_los = 2.0
        self.shadowing_sigma_dB_los = 4.1

        # model parameters for NLOS
        self.path_loss_exponent_nlos = 3.0
        self.shadowing_sigma_dB_nlos = 6.8


    def get_loss(self, *args, **kwargs) -> np.array:
        d = np.asarray(kwargs["distance"])
        f = np.asarray(kwargs["frequency"])
        p_los = kwargs["line_of_sight_prob"]
        std = kwargs["shadowing_std"]

        line_of_sight = np.random.sample(d.shape) <= p_los

        path_loss_exponent = np.ones(d.shape) * self._path_loss_exponent_nlos
        path_loss_exponent[line_of_sight] = self._path_loss_exponent_los

        shadowing_sigma_dB = np.ones(d.shape) * self._shadowing_sigma_dB_nlos
        shadowing_sigma_dB[line_of_sight] = self._shadowing_sigma_dB_los

        if std:
            shadowing = np.random.normal(0, std, d.shape)
        else:
            shadowing = 0

        loss = 20*np.log10(f) - 27.55 + 10*path_loss_exponent*np.log10(d)

        return loss + shadowing
