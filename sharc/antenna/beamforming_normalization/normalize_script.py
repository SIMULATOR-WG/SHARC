# -*- coding: utf-8 -*-
"""
Created on Mon Mar 12 13:26:45 2018

@author: Calil
"""

from sharc.support.named_tuples import AntennaPar
from sharc.antenna.beamforming_normalization.beamforming_normalizer import BeamformingNormalizer

#%% Setup

# General parameters
resolution = 2
tolerance = 1e-2

# Create object
norm = BeamformingNormalizer(resolution,tolerance)

# List of antenna parameters to which calculate the normalization factors.
# Order of parameters is: normalization element_pattern, element_max_g, 
# element_phi_deg_3db, element_theta_deg_3db, element_am, element_sla_v, 
# n_rows, n_columns, element_horiz_spacing, element_vert_spacing, downtilt_deg
param_list = [AntennaPar(True,"M2101",5,65,65,30,30,8,8,0.5,0.5,0),
              AntennaPar(True,"M2101",5,90,90,25,25,4,4,0.5,0.5,0),
              AntennaPar(True,"M2101",5,65,65,30,30,8,8,0.5,0.5,0),
              AntennaPar(True,"M2101",5,90,90,25,25,4,4,0.5,0.5,0)]
co_channel_list = [True, True, False, False]

#%% Normalize and save using parameters hash
for par, c_chan in zip(param_list,co_channel_list):
    file_name = str(hash((c_chan,par)))
    cf, err = norm.generate_correction_matrix(par,c_chan,file_name)