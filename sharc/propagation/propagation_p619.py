# -*- coding: utf-8 -*-
"""
Created on Mon Apr 17 15:35:00 2017

@author: andre barreto
"""

from sharc.propagation.propagation import Propagation

import numpy as np
import sys

class PropagationP619(Propagation):
    """
    Implements the earth-to-space channel model from ITU-R P.619
    Currently, only free-space is implemented
    """


    def get_loss(self, *args, **kwargs) -> np.array:
        d = np.asarray(kwargs["distance"])
        f = np.asarray(kwargs["frequency"])
        loss = 20*np.log10(d) + 20*np.log10(f) - 27.55
        return loss
