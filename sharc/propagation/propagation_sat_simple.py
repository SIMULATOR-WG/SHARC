# -*- coding: utf-8 -*-
"""
Created on Thu Feb 16 12:04:27 2017

@author: edgar
"""
import numpy as np
from multipledispatch import dispatch

from sharc.propagation.propagation import Propagation
from sharc.propagation.propagation_p619 import PropagationP619
from sharc.propagation.propagation_free_space import PropagationFreeSpace
from sharc.propagation.propagation_clutter_loss import PropagationClutterLoss
from sharc.propagation.propagation_building_entry_loss import PropagationBuildingEntryLoss
from sharc.support.enumerations import StationType
from sharc.station_manager import StationManager
from sharc.parameters.parameters import Parameters


class PropagationSatSimple(Propagation):
    """
    Implements the simplified satellite propagation model
    """
    # pylint: disable=function-redefined
    # pylint: disable=arguments-renamed

    def __init__(self, random_number_gen: np.random.RandomState, enable_clutter_loss=True):
        super().__init__(random_number_gen)
        self.enable_clutter_loss = enable_clutter_loss
        self.clutter = PropagationClutterLoss(random_number_gen)
        self.free_space = PropagationFreeSpace(random_number_gen)
        self.building_entry = PropagationBuildingEntryLoss(
            self.random_number_gen,
        )
        self.atmospheric_loss = 0.75

    @dispatch(Parameters, float, StationManager, StationManager, np.ndarray, np.ndarray)
    def get_loss(
        self,
        params: Parameters,
        frequency: float,
        station_a: StationManager,
        station_b: StationManager,
        station_a_gains=None,
        station_b_gains=None,
    ) -> np.array:
        """Wrapper function for the PropagationUMi calc_loss method
        Calculates the loss between station_a and station_b

        Parameters
        ----------
        station_a : StationManager
            StationManager container representing IMT UE station - Station_type.IMT_UE
        station_b : StationManager
            StationManager container representing IMT BS stattion
        params : Parameters
            Simulation parameters needed for the propagation class - Station_type.IMT_BS

        Returns
        -------
        np.array
            Return an array station_a.num_stations x station_b.num_stations with the path loss
            between each station
        """
        distance_3d = station_a.get_3d_distance_to(station_b)
        frequency = frequency * np.ones(distance_3d.shape)
        indoor_stations = np.tile(
            station_b.indoor, (station_a.num_stations, 1),
        )

        # Elevation angles seen from the station on Earth.
        elevation_angles = {}
        if station_a.is_space_station:
            elevation_angles["free_space"] = np.transpose(
                station_b.get_elevation(station_a),
            )
            elevation_angles["apparent"] = \
                PropagationP619.apparent_elevation_angle(
                    elevation_angles["free_space"],
                    station_a.height,)
        elif station_b.is_space_station:
            elevation_angles["free_space"] = station_a.get_elevation(station_b)
            elevation_angles["apparent"] = \
                PropagationP619.apparent_elevation_angle(
                    elevation_angles["free_space"],
                    station_b.height,)
        else:
            raise ValueError(
                "PropagationP619: At least one station must be an space station",
            )

        return self.get_loss(distance_3d, frequency, indoor_stations, elevation_angles)

    @dispatch(np.ndarray, np.ndarray, np.ndarray, dict)
    def get_loss(
        self,
        distance: np.array,
        frequency: np.array,
        indoor_stations: np.array,
        elevation: dict,
    ) -> np.array:
        """Calculates the clutter loss.

        Parameters
        ----------
        distance : np.array
            Distance between the stations
        frequency : np.array
            Array of frequenciews
        indoor_stations : np.array
            Bool array indicating if the terrestrial station is indoor or not.
        elevation : np.array
            Array with elevation angles w.r.t terrestrial station

        Returns
        -------
        np.array
            Array of clutter losses with the same shape as distance
        """

        free_space_loss = self.free_space.get_free_space_loss(
            distance=distance, frequency=frequency,
        )

        if self.enable_clutter_loss:
            clutter_loss = np.maximum(
                0, self.clutter.get_loss(
                    frequency=frequency,
                    distance=distance,
                    elevation=elevation["free_space"],
                    station_type=StationType.FSS_SS,
                ),
            )
        else:
            clutter_loss = 0

        building_loss = self.building_entry.get_loss(
            frequency, elevation["apparent"],
        ) * indoor_stations

        loss = free_space_loss + clutter_loss + building_loss + self.atmospheric_loss

        return loss
