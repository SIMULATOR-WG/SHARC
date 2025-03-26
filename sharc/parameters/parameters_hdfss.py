# -*- coding: utf-8 -*-
from dataclasses import dataclass

from sharc.parameters.parameters_base import ParametersBase


@dataclass
class ParametersHDFSS(ParametersBase):
    """
    Dataclass containing the HDFSS (High-Density Fixed Satellite System Propagation Model)
    propagation model parameters
    """
    # HDFSS position relative to building it is on. Possible values are
    # ROOFTOP and BUILDINGSIDE
    es_position: str = "ROOFTOP"
    # Enable shadowing loss
    shadow_enabled: bool = True
    # Enable building entry loss
    building_loss_enabled: bool = True
    # Enable interference from IMT stations at the same building as the HDFSS
    same_building_enabled: bool = False
    # Enable diffraction loss
    diffraction_enabled: bool = False
    # Building entry loss type applied between BSs and HDFSS ES. Options are:
    # P2109_RANDOM: random probability at P.2109 model, considering elevation
    # P2109_FIXED: fixed probability at P.2109 model, considering elevation.
    #              Probability must be specified in bs_building_entry_loss_prob.
    # FIXED_VALUE: fixed value per BS. Value must be specified in
    #              bs_building_entry_loss_value.
    bs_building_entry_loss_type: str = "P2109_FIXED"
    # Probability of building entry loss not exceeded if
    # bs_building_entry_loss_type = P2109_FIXED
    bs_building_entry_loss_prob: float = 0.75
    # Value in dB of building entry loss if
    # bs_building_entry_loss_type = FIXED_VALUE
    bs_building_entry_loss_value: float = 35
