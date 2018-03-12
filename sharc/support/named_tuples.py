# -*- coding: utf-8 -*-
"""
Created on Thu May 11 14:21:57 2017

@author: Calil
"""

from collections import namedtuple

AntennaPar = namedtuple("AntennaPar",
                        "normalization element_pattern element_max_g element_phi_deg_3db element_theta_deg_3db element_am element_sla_v n_rows n_columns element_horiz_spacing element_vert_spacing downtilt_deg")
