# -*- coding: utf-8 -*-
"""
Created on Thu Feb 16 12:03:12 2017

@author: edgar
"""

from abc import ABC, abstractmethod
import numpy as np
 
class Propagation(ABC):
    """
    Abstract base class for propagation models
    """
    
    @abstractmethod
    def get_loss(self, *args, **kwargs) -> np.array:
        pass
