# -*- coding: utf-8 -*-
"""
Created on Wed Feb 15 16:05:58 2017

@author: edgar
"""

class ParametersImt(object):

    static_base_stations = True
    num_clusters = 1
    intersite_distance = 1500       # [m]
    interfered_with = False
    frequency = 27250           # [MHz]

    bs_load_probability = 0.5
    bs_tx_power = 37          # [dBm]
    bs_height = 30

    ue_tx_power = 37          # [dBm]
    ue_height = 30