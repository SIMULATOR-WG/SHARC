# -*- coding: utf-8 -*-
"""
Created on Thu Feb 16 12:04:27 2017

@author: edgar
"""

from sharc.propagation.propagation import Propagation 

import numpy as np
 
class PropagationFreeSpace(Propagation):
    """
    Implements the Free Space propagation model
    """
    
    def __init__(self):
        super().__init__()
    
    def get_loss(self, *args, **kwargs) -> np.array:
        if "distance_2D" in kwargs:
            d = np.asarray(kwargs["distance_2D"])
        else:
            d = np.asarray(kwargs["distance_3D"])
        f = np.asarray(kwargs["frequency"])
        loss = 20*np.log10(d) + 20*np.log10(f) - 27.55
        return loss
