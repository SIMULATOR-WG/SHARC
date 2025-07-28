# -*- coding: utf-8 -*-
"""
Created on Tue Aug 28 13:57:48 2018

@author: Calil
"""

import numpy as np
import sys

from sharc.parameters.parameters_fss_es import ParametersFssEs
from sharc.propagation.propagation import Propagation
from sharc.station_manager import StationManager
from sharc.parameters.parameters import Parameters
from sharc.propagation.propagation_hdfss_roof_top import PropagationHDFSSRoofTop
from sharc.propagation.propagation_hdfss_building_side import PropagationHDFSSBuildingSide


class PropagationHDFSS(Propagation):
    """
    High-Density Fixed Satellite System Propagation Model
    This is a compoposition of HDFSS Rooftop and Indoor models using during simulation run-time.
    """

    def __init__(self, param: ParametersFssEs, rnd_num_gen: np.random.RandomState):
        super().__init__(rnd_num_gen)

        if param.es_position == "ROOFTOP":
            self.propagation = PropagationHDFSSRoofTop(param, rnd_num_gen)
        elif param.es_position == "BUILDINGSIDE":
            self.propagation = PropagationHDFSSBuildingSide(param, rnd_num_gen)
        else:
            sys.stderr.write(
                "ERROR\nInvalid es_position: " + param.es_position,
            )
            sys.exit(1)

    def get_loss(
        self,
        params: Parameters,
        frequency: float,
        station_a: StationManager,
        station_b: StationManager,
        station_a_gains=None,
        station_b_gains=None,
    ) -> np.array:
        """Wrapper function for the get_loss method to fit the Propagation ABC class interface
        Calculates the loss between station_a and station_b

        Parameters
        ----------
        params : Parameters
            Simulation parameters needed for the propagation class
        frequency: float
            Center frequency
        station_a : StationManager
            StationManager container representing the system station
        station_b : StationManager
            StationManager container representing the IMT station
        station_a_gains: np.ndarray defaults to None
            System antenna gains
        station_b_gains: np.ndarray defaults to None
            IMT antenna gains

        Returns
        -------
        np.array
            Return an array station_a.num_stations x station_b.num_stations with the path loss
            between each station
        """
        distance = station_a.get_3d_distance_to(station_b)  # P.452 expects Kms
        frequency_array = frequency * \
            np.ones(distance.shape)  # P.452 expects GHz
        elevation = station_b.get_elevation(station_a)

        return self.propagation.get_loss(
            distance_3D=distance,
            elevation=elevation,
            imt_sta_type=station_b.station_type,
            frequency=frequency_array,
            imt_x=station_b.x,
            imt_y=station_b.y,
            imt_z=station_b.height,
            es_x=station_a.x,
            es_y=station_a.y,
            es_z=station_a.height,
        )
