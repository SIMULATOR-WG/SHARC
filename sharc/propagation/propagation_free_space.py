# -*- coding: utf-8 -*-
"""
Created on Thu Feb 16 12:04:27 2017

@author: edgar
"""
import numpy as np
from multipledispatch import dispatch

from sharc.propagation.propagation import Propagation
from sharc.station_manager import StationManager
from sharc.parameters.parameters import Parameters


class PropagationFreeSpace(Propagation):
    """
    Implements the Free Space propagation model.

    Frequency in MHz and distance are in meters
    """

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
        """Wrapper for the calculation loss between station_a and station_b.

        Parameters
        ----------
        station_a : StationManager
            StationManager container representing station_a
        station_b : StationManager
            StationManager container representing station_a
        params : Parameters
            Simulation parameters needed for the propagation class.

        Returns
        -------
        np.array
            Return an array station_a.num_stations x station_b.num_stations with the path loss
            between each station
        """
        distance_3d = station_a.get_3d_distance_to(station_b)
        loss = self.get_free_space_loss(frequency=frequency, distance=distance_3d)

        return loss

    @dispatch(np.ndarray, np.ndarray)
    def get_loss(self, distance_3D: np.array, frequency: float) -> np.array:
        return self.get_free_space_loss(np.unique(frequency), distance_3D)

    def get_free_space_loss(self, frequency: float, distance: np.array) -> np.array:
        """Calculates the free-space loss for the given distance and frequency

        Parameters
        ----------
        distance : float
            3D distance array between stations
        frequency : float
            wave frequency
        Returns
        -------
        np.array
            returns the path loss array with shape distance.shape
        """
        loss = 20 * np.log10(distance) + 20 * np.log10(frequency) - 27.55

        return loss
