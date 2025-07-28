"""Generates the parameters for the MSS-SS to IMT with varying border distances campaign.
Parameters for the MSS_SS to IMT-DL interferece scenario in the 2300MHz band.
In this scenario the distance between the border of the MSS-SS footprint
and the IMT cluster is varied.
"""
import numpy as np
import yaml
import os

from sharc.parameters.parameters_base import tuple_constructor

yaml.SafeLoader.add_constructor('tag:yaml.org,2002:python/tuple', tuple_constructor)

local_dir = os.path.dirname(os.path.abspath(__file__))
parameter_file_name = os.path.join(local_dir, "../input/parameters_imt_ntn_to_imt_tn_sep_distance_template.yaml")

# load the base parameters from the yaml file
with open(parameter_file_name, 'r') as file:
    parameters_template = yaml.safe_load(file)

# Distance from topology boarders in meters
border_distances_array = np.array(
    [0, 10e3, 20e3, 30e3, 40e3, 50e3, 60e3, 70e3, 80e3, 90e3, 100e3])

for dist in border_distances_array:
    print(f'Generating parameters for distance {dist / 1e3} km')
    # Create a copy of the base parameters
    params = parameters_template.copy()

    ntn_footprint_left_edge = - 4 * params['mss_ss']['cell_radius']
    ntn_footprint_radius = 5 * params['mss_ss']['cell_radius'] * np.sin(np.pi / 3)
    macro_topology_radius = 4 * params['imt']['topology']['macrocell']['intersite_distance'] / 3
    params['mss_ss']['x'] = float(macro_topology_radius + ntn_footprint_radius + dist)

    # Set the right campaign prefix
    params['general']['output_dir_prefix'] = 'output_imt_ntn_to_imt_tn_co_channel_sep_' + str(dist / 1e3) + "_km"
    # Save the parameters to a new yaml file
    parameter_file_name = "../input/parameters_imt_ntn_to_imt_tn_co_channel_sep_" + str(dist / 1e3) + "_km.yaml"
    with open(os.path.join(local_dir, parameter_file_name), 'w') as file:
        yaml.dump(params, file, default_flow_style=False)

    print(f'Parameters saved to {parameter_file_name} file.')

    del params
