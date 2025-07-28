# -*- coding: utf-8 -*-
"""
Created on Sat Apr 15 16:29:36 2017

@author: Calil
"""

from sharc.support.named_tuples import AntennaPar
from numpy import load
import typing

from dataclasses import dataclass, field
from sharc.parameters.parameters_base import ParametersBase


@dataclass
class ParametersAntennaSubarrayImt(ParametersBase):
    """
    Parameters for subarray as defined in R23-WP5D-C-0413, Annex 4.2
    """
    # to use subarray, set this to true
    is_enabled: bool = False

    # Number of rows in subarray
    n_rows: int = 3

    # BS array element vertical spacing (d/lambda).
    element_vert_spacing: float = 0.5
    # element_vert_spacing: float = 0.5

    # notice that electrical tilt == -1 * downtilt
    # Sub array eletrical downtilt [deg]
    eletrical_downtilt: float = 3.0


@dataclass
class ParametersAntennaImt(ParametersBase):
    """
    Defines the antenna model and related parameters to be used in compatibility
    studies between IMT and other services in adjacent bands.
    """
    section_name: str = "imt_antenna"

    # Normalization application flags for base station (BS) and user equipment (UE).
    normalization: bool = False

    # Normalization files for BS and UE beamforming.
    normalization_file: str = "antenna/beamforming_normalization/norm.npz"

    # Radiation pattern of each antenna element.
    element_pattern: str = "M2101"

    # Minimum array gain for the beamforming antenna [dBi].
    minimum_array_gain: float = -200.0

    # beamforming angle limitation [deg].
    # PS: it isn't implemented for UEs
    # and current implementation doesn't make sense for UEs
    horizontal_beamsteering_range: tuple[float | int, float | int] = (-180., 180.)
    vertical_beamsteering_range: tuple[float | int, float | int] = (0., 180.)

    # Mechanical downtilt [degrees].
    # PS: downtilt doesn't make sense on UE's
    downtilt: float = 6.0

    # BS/UE maximum transmit/receive element gain [dBi].
    element_max_g: float = 5.0

    # BS/UE horizontal 3dB beamwidth of single element [degrees].
    element_phi_3db: float = 65.0

    # BS/UE vertical 3dB beamwidth of single element [degrees].
    element_theta_3db: float = 65.0

    # BS/UE number of rows and columns in antenna array.
    n_rows: int = 8
    n_columns: int = 8

    # BS/UE array element spacing (d/lambda).
    element_horiz_spacing: float = 0.5
    element_vert_spacing: float = 0.5

    # BS/UE front to back ratio and single element vertical sidelobe attenuation [dB].
    element_am: int = 30
    element_sla_v: int = 30

    # Multiplication factor k used to adjust the single-element pattern.
    multiplication_factor: int = 12

    adjacent_antenna_model: typing.Literal["BEAMFORMING", "SINGLE_ELEMENT"] = None

    subarray: ParametersAntennaSubarrayImt = field(default_factory=ParametersAntennaSubarrayImt)

    def load_subparameters(self, ctx: str, params: dict, quiet=True):
        """
        Loads the parameters when is placed as subparameter
        """
        super().load_subparameters(ctx, params, quiet)

    def set_external_parameters(self, *, adjacent_antenna_model: typing.Literal["BEAMFORMING", "SINGLE_ELEMENT"]):
        self.adjacent_antenna_model = adjacent_antenna_model

    def validate(self, ctx: str):
        # Additional sanity checks specific to antenna parameters can be implemented here

        # Sanity check for adjacent_antenna_model
        if self.adjacent_antenna_model not in ["SINGLE_ELEMENT", "BEAMFORMING"]:
            raise ValueError("adjacent_antenna_model must be 'SINGLE_ELEMENT'")

        # Sanity checks for normalization flags
        if not isinstance(self.normalization, bool):
            raise ValueError("normalization must be a boolean value")

        # Sanity checks for element patterns
        if self.element_pattern.upper() not in ["M2101", "F1336", "FIXED"]:
            raise ValueError(
                f"Invalid element_pattern value {self.element_pattern}",
            )
        if isinstance(self.horizontal_beamsteering_range, list):
            self.horizontal_beamsteering_range = tuple(self.horizontal_beamsteering_range)

        if not isinstance(self.horizontal_beamsteering_range, tuple):
            raise ValueError(
                f"Invalid {ctx}.horizontal_beamsteering_range={self.horizontal_beamsteering_range}\n"
                "It needs to be a tuple"
            )
        if len(self.horizontal_beamsteering_range) != 2\
            or not all(map(
                lambda x: isinstance(x, float) or isinstance(x, int), self.horizontal_beamsteering_range
            )):
            raise ValueError(
                f"Invalid {ctx}.horizontal_beamsteering_range={self.horizontal_beamsteering_range}\n"
                "It needs to contain two numbers delimiting the range of beamsteering in degrees"
            )
        if self.horizontal_beamsteering_range[0] > self.horizontal_beamsteering_range[1]:
            raise ValueError(
                f"Invalid {ctx}.horizontal_beamsteering_range={self.horizontal_beamsteering_range}\n"
                "The second value must be bigger than the first"
            )
        if not all(map(
                lambda x: x >= -180. and x <= 180., self.horizontal_beamsteering_range
            )):
            raise ValueError(
                f"Invalid {ctx}.horizontal_beamsteering_range={self.horizontal_beamsteering_range}\n"
                "Horizontal beamsteering limit angles must be in the range [-180, 180]"
            )

        if isinstance(self.vertical_beamsteering_range, list):
            self.vertical_beamsteering_range = tuple(self.vertical_beamsteering_range)
        if not isinstance(self.vertical_beamsteering_range, tuple):
            raise ValueError(
                f"Invalid {ctx}.vertical_beamsteering_range={self.vertical_beamsteering_range}\n"
                "It needs to be a tuple"
            )
        if len(self.vertical_beamsteering_range) != 2\
            or not all(map(
                lambda x: isinstance(x, float) or isinstance(x, int), self.vertical_beamsteering_range
            )):
            raise ValueError(
                f"Invalid {ctx}.vertical_beamsteering_range={self.vertical_beamsteering_range}\n"
                "It needs to contain two numbers delimiting the range of beamsteering in degrees"
            )
        if self.vertical_beamsteering_range[0] > self.vertical_beamsteering_range[1]:
            raise ValueError(
                f"Invalid {ctx}.vertical_beamsteering_range={self.vertical_beamsteering_range}\n"
                "The second value must be bigger than the first"
            )
        if not all(map(
                lambda x: x >= 0. and x <= 180., self.vertical_beamsteering_range
            )):
            raise ValueError(
                f"Invalid {ctx}.vertical_beamsteering_range={self.vertical_beamsteering_range}\n"
                "vertical beamsteering limit angles must be in the range [0, 180]"
            )

    def get_antenna_parameters(self) -> AntennaPar:
        if self.normalization:
            # Load data, save it in dict and close it
            data = load(self.normalization_file)
            data_dict = {key: data[key] for key in data}
            self.normalization_data = data_dict
            data.close()
        else:
            self.normalization_data = None
        tpl = AntennaPar(
            self.adjacent_antenna_model,
            self.normalization,
            self.normalization_data,
            self.element_pattern,
            self.element_max_g,
            self.element_phi_3db,
            self.element_theta_3db,
            self.element_am,
            self.element_sla_v,
            self.n_rows,
            self.n_columns,
            self.element_horiz_spacing,
            self.element_vert_spacing,
            self.multiplication_factor,
            self.minimum_array_gain,
            self.downtilt,
        )

        return tpl
