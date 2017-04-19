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
    def __init__(self, line_of_sight_prob: float = 0 ):
        self.__nlos_loss = 20
        self.__atmospheric_loss = 1
        self.__polarization_loss = 3
        self.__line_of_sight_prob = line_of_sight_prob

    def get_loss(self, *args, **kwargs) -> np.array:
        d = np.asarray(kwargs["distance"])
        f = np.asarray(kwargs["frequency"])
        free_space_loss = 20*np.log10(d) + 20*np.log10(f) - 27.55
        line_of_sight = np.random.sample(d.shape) <= self.__line_of_sight_prob
        nlos_loss = self.__nlos_loss * ( ~line_of_sight )
        loss = (free_space_loss + nlos_loss + self.__polarization_loss
                + self.__atmospheric_loss)
        return loss
