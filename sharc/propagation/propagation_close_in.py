# -*- coding: utf-8 -*-
"""
Created on Wed Apr 12 16:59:27 2017

@author: andre barreto
"""

from sharc.propagation.propagation import Propagation

import numpy as np
import sys

class PropagationCloseIn(Propagation):
    """
    Implements the close-in free-space reference distance propagation model
    Initially only the macro-cell topology is supported
    """

    def __init__(self, topology: str = "MACROCELL", line_of_sight_prob: float = 0,
                 enable_shadowing: bool = True ):

        self._topology = topology

        self._line_of_sight_prob = line_of_sight_prob

        # model parameters for LOS
        self._path_loss_exponent_los = 2.0

        # model parameters for NLOS
        self._path_loss_exponent_nlos = 3.0

        if not (topology == "MACROCELL" or topology == "SINGLE_BS"):
            sys.stderr.write("error: class PropagationCloseIn: " + topology +
                             " topology not supported")
            sys.exit(1)

        self._enable_shadowing = enable_shadowing
        self._shadowing_sigma_dB_los = 4.1
        self._shadowing_sigma_dB_nlos = 6.8


    def get_loss(self, *args, **kwargs) -> np.array:
        d = np.asarray(kwargs["distance"])
        f = np.asarray(kwargs["frequency"])

        number_of_links = d.shape

        line_of_sight = np.random.sample(number_of_links) <= self._line_of_sight_prob

        path_loss_exponent = np.ones(number_of_links) * self._path_loss_exponent_nlos
        path_loss_exponent[line_of_sight] = self._path_loss_exponent_los

        shadowing_sigma_dB = np.ones(number_of_links) * self._shadowing_sigma_dB_nlos
        shadowing_sigma_dB[line_of_sight] = self._shadowing_sigma_dB_los

        if self._enable_shadowing:
            shadowing = np.random.normal(0, shadowing_sigma_dB, d.shape )
        else:
            shadowing = 0

        loss = ( 20*np.log10(f) - 27.55 + 10*path_loss_exponent*np.log10( d )
                 + shadowing )

        return loss
