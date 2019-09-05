# -*- coding: utf-8 -*-
"""
Created on Thu May 11 14:21:57 2017

@author: Calil
"""

from collections import namedtuple

AntennaPar = namedtuple("AntennaPar",
                        "adjacent_antenna_model normalization normalization_data element_pattern element_max_g element_phi_3db element_theta_3db element_am element_sla_v n_rows n_columns element_horiz_spacing element_vert_spacing multiplication_factor minimum_array_gain downtilt")
