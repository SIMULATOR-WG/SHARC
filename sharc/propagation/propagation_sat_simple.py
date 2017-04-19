# -*- coding: utf-8 -*-
"""
Created on Thu Feb 16 12:04:27 2017

@author: edgar
"""

from sharc.propagation.propagation import Propagation

import numpy as np

class PropagationSatSimple(Propagation):
    """
    Implements the simplified satellite propagation model
    """
    def __init__(self, line_of_sight: bool = False ):
        if not line_of_sight:
            self.__nlos_loss = 20
        else:
            self.__nlos_loss = 0

        self.__atmospheric_loss = 1
        self.__polarization_loss = 3

    def get_loss(self, *args, **kwargs) -> np.array:
        d = np.asarray(kwargs["distance"])
        f = np.asarray(kwargs["frequency"])
        free_space_loss = 20*np.log10(d) + 20*np.log10(f) - 27.55
        loss = (free_space_loss + self.__nlos_loss + self.__polarization_loss
                + self.__atmospheric_loss)
        return loss
