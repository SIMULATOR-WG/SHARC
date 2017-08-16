# -*- coding: utf-8 -*-
"""
Created on Thu Jul  6 16:03:24 2017

@author: edgar
"""

import sys

from sharc.propagation.propagation import Propagation
from sharc.propagation.propagation_free_space import PropagationFreeSpace
from sharc.propagation.propagation_close_in import PropagationCloseIn
from sharc.propagation.propagation_p619 import PropagationP619
from sharc.propagation.propagation_sat_simple import PropagationSatSimple
from sharc.propagation.propagation_ter_simple import PropagationTerSimple
from sharc.propagation.propagation_uma import PropagationUMa
from sharc.propagation.propagation_umi import PropagationUMi
from sharc.propagation.propagation_abg import PropagationABG

class PropagationFactory(object):
    
    @staticmethod
    def createPropagation(channel_model: str) -> Propagation:
        if channel_model == "FSPL":
            return PropagationFreeSpace()
        elif channel_model == "ABG":
            return PropagationABG()
        elif channel_model == "UMa":
            return PropagationUMa()
        elif channel_model == "UMi":
            return PropagationUMi()            
        elif channel_model == "CI":
            return PropagationCloseIn()
        elif channel_model == "SatelliteSimple":
            return PropagationSatSimple()
        elif channel_model == "TerrestrialSimple":
            return PropagationTerSimple()
        elif channel_model == "P619":
            return PropagationP619()
        else:
            sys.stderr.write("ERROR\nInvalid channel_model: " + channel_model)
            sys.exit(1)
