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
    
    def __init__(self):
        super().__init__()
        self.nlos_loss = 20
        self.atmospheric_loss = 1
        self.polarization_loss = 3

        
    def get_loss(self, *args, **kwargs) -> np.array:
        d = np.asarray(kwargs["distance_3D"])
        f = np.asarray(kwargs["frequency"])
        p_los = kwargs["line_of_sight_prob"]
        
        free_space_loss = 20*np.log10(d) + 20*np.log10(f) - 27.55
        line_of_sight = np.random.sample(d.shape) <= p_los
        nlos_loss = self.nlos_loss * ( ~line_of_sight )
        loss = (free_space_loss + nlos_loss + self.polarization_loss
                + self.atmospheric_loss)
        
        return loss
