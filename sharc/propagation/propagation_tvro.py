# -*- coding: utf-8 -*-
"""
Created on Mon Jun  5 16:56:13 2017

@author: edgar
"""

from sharc.propagation.propagation import Propagation
from sharc.propagation.propagation_free_space import PropagationFreeSpace

import numpy as np
import matplotlib.pyplot as plt
from cycler import cycler

class PropagationTvro(Propagation):
    """
    Implements the propagation model used in paper
    Fernandes, Linhares, "Coexistence conditions of LTE-advanced at 3400-3600MHz with TVRO
                          at 3625-4200 MHz in Brazil", Wireless Networks, 2017
    TODO: calculate the effective environment height for the generic case
    """

    def __init__(self, random_number_gen: np.random.RandomState):
        super().__init__(random_number_gen)

        self.d_k = 0.2 #km
        self.std = 6.
        self.h_a = 20

        self.free_space_path_loss = PropagationFreeSpace(random_number_gen)

    def get_loss(self, *args, **kwargs) -> np.array:
        """
        Calculates path loss

        Parameters
        ----------
            distance_3D (np.array) : 3D distances between stations
            distance_2D (np.array) : 2D distances between stations
            frequency (np.array) : center frequencie [MHz]
            bs_height (np.array) : base station antenna heights
        Returns
        -------
            array with path loss values with dimensions of distance_2D

        """
        d_3D = kwargs["distance_3D"]
        f_MHz = kwargs["frequency"]
        height = kwargs["es_params"].height
        f_GHz = f_MHz / 1000
        number_of_sectors = kwargs.pop("number_of_sectors",1)

        free_space_path_loss = self.free_space_path_loss.get_loss(distance_3D=d_3D, frequency=f_MHz)
        shadowing = self.random_number_gen.randn(d_3D.size)

        f_fc = .25 + .375*(1 + np.tanh(7.5*(f_GHz-.5)))
        clutter_loss = 10.25 * f_fc * np.exp(-self.d_k) * \
                       (1 - np.tanh(6*(height/self.h_a - .625))) - .33

        loss = free_space_path_loss + shadowing

        indices = (d_3D >= 40) & (d_3D < 10 * self.d_k * 1000)
        loss[indices] = loss[indices] + (d_3D[indices] - 40)/(10 * 1000 * self.d_k - 40) * clutter_loss[indices]

        indices = (d_3D >= 10 * self.d_k * 1000)
        loss[indices] = loss[indices] + clutter_loss[indices]

        if number_of_sectors > 1:
            loss = np.repeat(loss, number_of_sectors, 1)

        return loss

