# -*- coding: utf-8 -*-
"""
Created on Thu Feb 16 12:03:12 2017

@author: edgar
"""

from abc import ABC, abstractmethod
import numpy as np

from sharc.station_manager import StationManager
from sharc.parameters.parameters import Parameters


class Propagation(ABC):
    """
    Abstract base class for propagation models
    """

    def __init__(self, random_number_gen: np.random.RandomState):
        self.random_number_gen = random_number_gen
        # Inicates whether this propagation model is for links between earth and space
        self.is_earth_space_model = False

    @abstractmethod
    def get_loss(
        self,
        params: Parameters,
        frequency: float,
        station_a: StationManager,
        station_b: StationManager,
        station_a_gains=None,
        station_b_gains=None,
    ) -> np.array:
        """Calculates the loss between station_a and station_b

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
