"""Generates scenarios based on main parameters
"""
import numpy as np
import yaml
import os
from copy import deepcopy

from sharc.parameters.parameters_base import tuple_constructor

yaml.SafeLoader.add_constructor('tag:yaml.org,2002:python/tuple', tuple_constructor)

local_dir = os.path.dirname(os.path.abspath(__file__))
input_dir = os.path.join(local_dir, "../input")

ul_parameter_file_name = os.path.join(local_dir, "./base_input.yaml")

# load the base parameters from the yaml file
with open(ul_parameter_file_name, 'r') as file:
    ul_parameters = yaml.safe_load(file)

dl_parameters = deepcopy(ul_parameters)
dl_parameters['general']['output_dir'] = ul_parameters['general']['output_dir'].replace("_ul", "_dl")
dl_parameters['general']['output_dir_prefix'] = ul_parameters['general']['output_dir_prefix'].replace("_ul", "_dl")
dl_parameters['general']['imt_link'] = "DOWNLINK"

country_border = 4 * ul_parameters["mss_d2d"]["cell_radius"] / 1e3
print("country_border", country_border)

# doesn't matter from which, both will give same result
output_dir_pattern = ul_parameters['general']['output_dir'].replace("_base_ul", "_<specific>")
output_prefix_pattern = ul_parameters['general']['output_dir_prefix'].replace("_base_ul", "_<specific>")

for dist in [
    0,
    country_border,
    country_border + 111 / 2,
    country_border + 111,
    country_border + 3 * 111 / 2,
    country_border + 2 * 111
]:
    ul_parameters["mss_d2d"]["sat_is_active_if"]["lat_long_inside_country"]["margin_from_border"] = dist
    specific = f"{dist}km_base_ul"
    ul_parameters['general']['output_dir_prefix'] = output_prefix_pattern.replace("<specific>", specific)

    dl_parameters["mss_d2d"]["sat_is_active_if"]["lat_long_inside_country"]["margin_from_border"] = dist
    specific = f"{dist}km_base_dl"
    dl_parameters['general']['output_dir_prefix'] = output_prefix_pattern.replace("<specific>", specific)

    ul_parameter_file_name = os.path.join(input_dir, f"./parameters_mss_d2d_to_imt_cross_border_{dist}km_base_ul.yaml")
    dl_parameter_file_name = os.path.join(input_dir, f"./parameters_mss_d2d_to_imt_cross_border_{dist}km_base_dl.yaml")

    with open(
        dl_parameter_file_name,
        'w'
    ) as file:
        yaml.dump(dl_parameters, file, default_flow_style=False)
    with open(
        ul_parameter_file_name,
        'w'
    ) as file:
        yaml.dump(ul_parameters, file, default_flow_style=False)

    for link in ["ul", "dl"]:
        if link == "ul":
            parameters = deepcopy(ul_parameters)
        if link == "dl":
            parameters = deepcopy(dl_parameters)

        parameters['mss_d2d']['num_sectors'] = 19
        # 1 out of 19 beams are active
        parameters['mss_d2d']['beams_load_factor'] = 0.05263157894

        specific = f"{dist}km_activate_random_beam_5p_{link}"
        parameters['general']['output_dir_prefix'] = output_prefix_pattern.replace("<specific>", specific)

        with open(
            os.path.join(input_dir, f"./parameters_mss_d2d_to_imt_cross_border_{specific}.yaml"),
            'w'
        ) as file:
            yaml.dump(parameters, file, default_flow_style=False)

        parameters['mss_d2d']['num_sectors'] = 19
        parameters['mss_d2d']['beams_load_factor'] = 0.3

        specific = f"{dist}km_activate_random_beam_30p_{link}"
        parameters['general']['output_dir_prefix'] = output_prefix_pattern.replace("<specific>", specific)

        with open(
            os.path.join(input_dir, f"./parameters_mss_d2d_to_imt_cross_border_{specific}.yaml"),
            'w'
        ) as file:
            yaml.dump(parameters, file, default_flow_style=False)

        parameters['mss_d2d']['num_sectors'] = 1
        parameters['mss_d2d']['beams_load_factor'] = 1

        parameters['mss_d2d']['center_beam_positioning'] = {}

        parameters['mss_d2d']['center_beam_positioning']['type'] = "ANGLE_AND_DISTANCE_FROM_SUBSATELLITE"

        # for uniform area distribution
        parameters['mss_d2d']['center_beam_positioning']['angle_from_subsatellite_phi'] = {
            'type': "~U(MIN,MAX)",
            'distribution': {
                'min': -180.,
                'max': 180.,
            }
        }
        parameters['mss_d2d']['center_beam_positioning']['distance_from_subsatellite'] = {
            'type': "~SQRT(U(0,1))*MAX",
            'distribution': {
                'min': 0,
                'max': parameters['mss_d2d']["cell_radius"] * 4,
            }
        }

        specific = f"{dist}km_random_pointing_1beam_{link}"
        parameters['general']['output_dir_prefix'] = output_prefix_pattern.replace("<specific>", specific)

        with open(
            os.path.join(input_dir, f"./parameters_mss_d2d_to_imt_cross_border_{specific}.yaml"),
            'w'
        ) as file:
            yaml.dump(parameters, file, default_flow_style=False)
