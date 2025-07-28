"""Generates the parameters for the MSS D2D to IMT with varying latitude campaign.
"""
import numpy as np
import yaml
import os

from sharc.parameters.parameters_base import tuple_constructor

yaml.SafeLoader.add_constructor('tag:yaml.org,2002:python/tuple', tuple_constructor)

local_dir = os.path.dirname(os.path.abspath(__file__))
parameter_file_name = os.path.join(local_dir, "../input/parameters_mss_d2d_to_imt_lat_variation_template.yaml")

# load the base parameters from the yaml file
with open(parameter_file_name, 'r') as file:
    parameters = yaml.safe_load(file)

# The scenario is to move the Earth station postion throught the Meridian from 0 to 60 degrees latitude.
parameters['imt']['topology']['central_longitude'] = 0.0
lats = np.arange(0, 70, 10)
for direction in [('dl', 'DOWNLINK'), ('ul', 'UPLINK')]:
    for lat in lats:
        # Create a copy of the base parameters
        params = parameters.copy()

        # set the link direction
        params['general']['imt_link'] = direction[1]

        # Update the parameters with the new latitude
        params['imt']['topology']['central_latitude'] = float(lat)
        params['imt']['topology']['central_longitude'] = 0.0

        # Set the right campaign prefix
        params['general']['output_dir_prefix'] = 'output_mss_d2d_to_imt_lat_' + direction[0] + '_' + str(lat) + "_deg"
        # Save the parameters to a new yaml file
        parameter_file_name = "../input/parameters_mss_d2d_to_imt_lat_" + direction[0] + '_' + str(lat) + "_deg.yaml"
        parameter_file_name = os.path.join(local_dir, parameter_file_name)
        with open(parameter_file_name, 'w') as file:
            yaml.dump(params, file, default_flow_style=False)

        del params
        print(f'Generated parameters for latitude {lat} degrees.')
