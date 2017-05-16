# -*- coding: utf-8 -*-
"""
Created on Thu Feb 16 12:03:12 2017

@author: edgar
"""

from abc import ABCMeta, abstractmethod
import numpy as np
 
class Propagation(object):
    """
    Abstract base class for propagation models
    """
    
    __metaclass__ = ABCMeta
    
    @abstractmethod
    def get_loss(self, *args, **kwargs) -> np.array:
        pass
    
    @abstractmethod
    def get_loss_Ag(self, *args, **kwargs) -> np.array:
        pass