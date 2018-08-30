# -*- coding: utf-8 -*-
"""
Created on Tue Aug 28 13:57:48 2018

@author: Calil
"""

import numpy as np
import sys

from sharc.parameters.parameters_fss_es import ParametersFssEs
from sharc.propagation.propagation import Propagation
from sharc.propagation.propagation_hdfss_roof_top import PropagationHDFSSRoofTop
from sharc.propagation.propagation_hdfss_biulding_side import PropagationHDFSSBuildingSide

class PropagationHDFSS(Propagation):
    """
    
    """
    def __init__(self, param: ParametersFssEs, rnd_num_gen: np.random.RandomState):
        """
        
        """
        super().__init__(rnd_num_gen)
        
        if param.es_position == "ROOFTOP":
            self.propagation = PropagationHDFSSRoofTop(param,rnd_num_gen)
        elif param.es_position == "BUILDINGSIDE":
            self.propagation = PropagationHDFSSBuildingSide(param,rnd_num_gen)
        else:
            sys.stderr.write("ERROR\nInvalid es_position: " + param.es_position)
            sys.exit(1)
        
    def get_loss(self, *args, **kwargs) -> np.array:
        """
        
        """
        if "distance_3D" in kwargs:
            d = kwargs["distance_3D"]
        else:
            d = kwargs["distance_2D"]
            
        ele = np.transpose(kwargs["elevation"])
        sta_type = kwargs["imt_sta_type"]
        f = kwargs["frequency"]
        num_sec = kwargs.pop("number_of_sectors",1)
        
        i_x = kwargs['imt_x']
        i_y = kwargs['imt_y']
        i_z = kwargs['imt_z']
        e_x = kwargs["es_x"]
        e_y = kwargs["es_y"]
        e_z = kwargs["es_z"]
        
        return self.propatation(distance_3D = d,
                                elevation = ele,
                                imt_sta_type = sta_type,
                                frequency = f,
                                number_of_sectors = num_sec,
                                imt_x = i_x,
                                imt_y = i_y,
                                imt_z = i_z,
                                es_x = e_x,
                                es_y = e_y,
                                es_z = e_z)
        