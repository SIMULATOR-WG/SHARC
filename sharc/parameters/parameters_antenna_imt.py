# -*- coding: utf-8 -*-
"""
Created on Sat Apr 15 16:29:36 2017

@author: Calil
"""

import sys

from sharc.support.named_tuples import AntennaPar
from sharc.support.enumerations import StationType
from numpy import load

class ParametersAntennaImt(object):
    """
    Defines parameters for antenna array.
    """

    def __init__(self):
        pass


    ###########################################################################
    # Named tuples which contain antenna types

    def get_antenna_parameters(self, sta_type: StationType)-> AntennaPar:
        if sta_type is StationType.IMT_BS:
            if self.bs_normalization:
                # Load data, save it in dict and close it
                data = load(self.bs_normalization_file)
                data_dict = {key:data[key] for key in data}
                self.bs_normalization_data = data_dict
                data.close()
            else:
                self.bs_normalization_data = None
            tpl = AntennaPar(self.adjacent_antenna_model,
                             self.bs_normalization,
                             self.bs_normalization_data,
                             self.bs_element_pattern,
                             self.bs_element_max_g,
                             self.bs_element_phi_3db,
                             self.bs_element_theta_3db,
                             self.bs_element_am,
                             self.bs_element_sla_v,
                             self.bs_n_rows,
                             self.bs_n_columns,
                             self.bs_element_horiz_spacing,
                             self.bs_element_vert_spacing,
                             self.bs_multiplication_factor,
                             self.bs_minimum_array_gain,
                             self.bs_downtilt)
        elif sta_type is StationType.IMT_UE:
            if self.ue_normalization:
                # Load data, save it in dict and close it
                data = load(self.ue_normalization_file)
                data_dict = {key:data[key] for key in data}
                self.ue_normalization_data = data_dict
                data.close()
            else:
                self.ue_normalization_data = None
            tpl = AntennaPar(self.adjacent_antenna_model,
                             self.ue_normalization,
                             self.ue_normalization_data,
                             self.ue_element_pattern,
                             self.ue_element_max_g,
                             self.ue_element_phi_3db,
                             self.ue_element_theta_3db,
                             self.ue_element_am,
                             self.ue_element_sla_v,
                             self.ue_n_rows,
                             self.ue_n_columns,
                             self.ue_element_horiz_spacing,
                             self.ue_element_vert_spacing,
                             self.ue_multiplication_factor,
                             self.ue_minimum_array_gain,
                             0)
        else:
            sys.stderr.write("ERROR\nInvalid station type: " + sta_type)
            sys.exit(1)

        return tpl
