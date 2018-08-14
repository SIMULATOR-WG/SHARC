# -*- coding: utf-8 -*-
"""
Created on Tue Jul  3 17:13:33 2018

@author: Calil
"""

from sharc.propagation.propagation import Propagation

import numpy as np

class PropagationP1411(Propagation):
    """
    Implements the propagation model described in ITU-R P.1411-9, section 4.2
    
    Frequency in MHz and distance in meters!
    """
    def __init__(self, random_number_gen: np.random.RandomState):
        super().__init__(random_number_gen)
        
        self.los_alpha = 2.29 
        self.los_beta = 28.6 
        self.los_gamma = 1.96 
        self.los_sigma = 3.48
        
        self.nlos_alpha = 4.39
        self.nlos_beta = -6.27
        self.nlos_gamma = 2.30
        self.nlos_sigma = 6.89
    
    def get_loss(self, *args, **kwargs) -> np.array:
        if "distance_3D" in kwargs:
            d = kwargs["distance_3D"]
        else:
            d = kwargs["distance_2D"]

        f = kwargs["frequency"]/1e3
        los = kwargs.pop("los",True)
        shadow = kwargs.pop("shadow",True)
        number_of_sectors = kwargs.pop("number_of_sectors",1)
        
        if los:
            alpha = self.los_alpha 
            beta = self.los_beta 
            gamma = self.los_gamma 
            sigma = self.los_sigma
        else:
            alpha = self.nlos_alpha
            beta = self.nlos_beta
            gamma = self.nlos_gamma
            sigma = self.nlos_sigma
            
        if shadow:
            shadow_loss = self.random_number_gen.normal(0.0,sigma,d.shape)
        else:
            shadow_loss = 0.0

        loss = 10*alpha*np.log10(d) + 10*gamma*np.log10(f) + beta + shadow_loss

        if number_of_sectors > 1:
            loss = np.repeat(loss, number_of_sectors, 1)

        return loss
