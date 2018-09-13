# -*- coding: utf-8 -*-
"""
Created on Mon Mar 12 13:26:45 2018

@author: Calil

This script generates the correction factors for the IMT Beamforming Antennas,
both array and single element, and saves them in files with the given names.
This script must be ran with the appropriate parameters prior to using any
normalization in the SHARC simulator, since the simulator merely reads the
correction factor values from the saved files.
For the co-channel scenario (antenna array) the correction factor is a 2-D
array with the lines representing the azimuth and the columns representing the
elevation of the beam direction.
For the adjacent channel scenario (single element) the correction factor is a
float.

Variables:
    resolution (float): resolution of the azimuth and elevation angles in the
        antenna array correction factor matrix [deg]. This defines the number
        of beam pointing directions to which the correction factor is
        calculated.
    tolerance (float): absolute tolerance of the correction factor integral, in
        linear scale.
    norm (BeamformingNormalizer): object that calculates the normalization.
    param_list (list): list of antenna parameters to which calculate the
        correction factors. New parameters are added as:
            AntennaPar(normalization,
                       norm_data,
                       element_pattern,
                       element_max_g,
                       element_phi_deg_3db,
                       element_theta_deg_3db,
                       element_am,
                       element_sla_v,
                       n_rows,
                       n_columns,
                       element_horiz_spacing,
                       element_vert_spacing,
                       downtilt_deg)
            normalization parameter must be set to False, otherwise script will
            try to normalize an already normalized antenna.
    file_names (list): list of file names to which save the normalization data.
        Files are paired with AntennaPar objects in param_list, so that the
        normalization data of the first element of param_list is saved in a
        file with the name specified in the first element of file_names and so
        on.

Data is saved in an .npz file in a dict like data structure with the
following keys:
    resolution (float): antenna array correction factor matrix angle resolution
        [deg]
    phi_range (tuple): range of beam pointing azimuth angle values [deg]
    theta_range (tuple): range of beam pointing elevation angle values [deg]
    correction_factor_co_channel (2D np.array): correction factor [dB]for the
        co-channel scenario (antenna array) for each of the phi theta pairs in
        phi_range and theta_range. Phi is associated with the lines and Theta
        is associated with the columns of the array.
    error_co_channel (2D np.array of tuples): lower and upper bounds of
        calculated correction factors [dB], considering integral error
    correction_factor_adj_channel (float):correction factor for single antenna
        element
    error_adj_channel (tuple): lower and upper bounds [dB] of single antenna
        element correction factor
    parameters (AntennaPar): antenna parameters used in the normalization
"""
from sharc.support.named_tuples import AntennaPar
from sharc.antenna.beamforming_normalization.beamforming_normalizer import BeamformingNormalizer

###############################################################################
## List of antenna parameters to which calculate the normalization factors.
param_list = [AntennaPar(False,None,"M2101",5,90,90,25,25,8,16,0.5,0.5,0)]
file_names = ['bs_indoor_norm.npz']
###############################################################################
## Setup
# General parameters
resolution = 2
tolerance = 1e-2

# Create object
norm = BeamformingNormalizer(resolution,tolerance)
###############################################################################
## Normalize and save
for par, file in zip(param_list,file_names):
    s = 'Generating ' + file
    print(s)

    norm.generate_correction_matrix(par,file)
