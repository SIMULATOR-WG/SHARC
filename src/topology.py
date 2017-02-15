# -*- coding: utf-8 -*-
"""
Created on Tue Feb 14 12:48:58 2017

@author: edgar
"""

from abc import ABCMeta, abstractmethod
import numpy as np
import typing
 
class Topology(object):
    
    __metaclass__ = ABCMeta
    
    @abstractmethod
    def get_coordinates(self) -> typing.Tuple[np.array, np.array]:
        pass
    