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
    Initially only the macro-cell, non LOS is implemented
    """

    def __init__(self, topology: str = "MACROCELL", line_of_sight: bool = False,
                       enable_shadowing: bool = True ):

        self.topology = topology

        if topology == "MACROCELL" or topology == "SINGLE_BS":
            if line_of_sight:
                self._path_loss_exponent = 2.0
                self._shadowing_sigma_dB = 4.1
            else:
                self._path_loss_exponent = 3.0
                self._shadowing_sigma_dB = 6.8
        else:
            sys.stderr.write("error: class PropagationCloseIn: " + topology +
                             " topology not supported")
            sys.exit(1)

        if not enable_shadowing:
            self._shadowing_sigma_dB = 0

    def get_loss(self, *args, **kwargs) -> np.array:
        d = np.asarray(kwargs["distance"])
        f = np.asarray(kwargs["frequency"])

        if self._shadowing_sigma_dB != 0:
            shadowing = np.random.normal(0, self._shadowing_sigma_dB, d.shape )
        else:
            shadowing = 0

        loss = ( 20*np.log10(f) - 27.55
                 + 10 * self._path_loss_exponent * np.log10( d )
                 + shadowing )

        return loss
