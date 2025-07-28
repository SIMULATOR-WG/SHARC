# -*- coding: utf-8 -*-
"""
Created on Wed Aug 16 13:42:19 2017

@author: edgar
"""
import numpy as np
from multipledispatch import dispatch

from sharc.propagation.propagation import Propagation
from sharc.station_manager import StationManager
from sharc.parameters.parameters import Parameters
from sharc.propagation.propagation_free_space import PropagationFreeSpace
from sharc.propagation.propagation_clutter_loss import PropagationClutterLoss
from sharc.support.enumerations import StationType


class PropagationTerSimple(Propagation):
    # pylint: disable=function-redefined
    # pylint: disable=arguments-renamed
    """
    Implements the simplified terrestrial propagation model, which is the
    basic free space and additional clutter losses
    """

    def __init__(self, random_number_gen: np.random.RandomState):
        super().__init__(random_number_gen)
        self.clutter = PropagationClutterLoss(np.random.RandomState(101))
        self.free_space = PropagationFreeSpace(np.random.RandomState(101))
        self.building_loss = 20

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
            Not used
        station_b_gains: np.ndarray defaults to None
            Not used

        Returns
        -------
        np.array
            Return an array station_a.num_stations x station_b.num_stations with the path loss
            between each station
        """
        distance = station_a.get_3d_distance_to(station_b)
        frequency_array = frequency * np.ones(distance.shape)
        indoor_stations = np.tile(
            station_b.indoor, (station_a.num_stations, 1),
        )

        return self.get_loss(distance, frequency_array, indoor_stations, -1.0)

    # pylint: disable=arguments-differ
    @dispatch(np.ndarray, np.ndarray, np.ndarray, float)
    def get_loss(
        self, distance: np.ndarray, frequency: np.ndarray, indoor_stations: np.ndarray,
        loc_percentage: float,
    ) -> np.array:
        """Calculates loss with a simple terrestrial model:
        * Building loss is set statically in class construction
        * Clutter loss is calculated using P.2108 standard
        * Free-space loss

        Parameters
        ----------
        distance : np.ndarray
            distances array between stations
        frequency : np.ndarray
            frequency
        indoor_stations: np.ndarray
            array of bool indicating if station n is indoor
        loc_percentage : str, optional
            Percentage locations range [0, 1[. If a negative number is given
            a random percentage is used.

        Returns
        -------
        np.array
            array of losses with distance dimentions
        """
        free_space_loss = self.free_space.get_free_space_loss(
            frequency, distance,
        )
        if loc_percentage < 0:
            loc_percentage = "RANDOM"

        clutter_loss = self.clutter.get_loss(
            frequency=frequency,
            distance=distance,
            loc_percentage=loc_percentage,
            station_type=StationType.FSS_ES,
        )

        building_loss = self.building_loss * indoor_stations

        loss = free_space_loss + building_loss + clutter_loss

        return loss


if __name__ == '__main__':

    import matplotlib.pyplot as plt

    ###########################################################################
    # Print path loss for TerrestrialSimple and Free Space

    d = np.linspace(10, 10000, num=10000)
    freq = 27000 * np.ones(d.shape)
    indoor_stations = np.zeros(d.shape, dtype=bool)
    loc_percentage = 0.5

    free_space = PropagationFreeSpace(np.random.RandomState(101))
    ter_simple = PropagationTerSimple(np.random.RandomState(101))

    loss_ter = ter_simple.get_loss(d, freq, indoor_stations, loc_percentage)

    loss_fs = free_space.get_free_space_loss(freq, d)

    fig = plt.figure(figsize=(8, 6), facecolor='w', edgecolor='k')

    plt.semilogx(np.squeeze(d), np.squeeze(loss_fs), label="free space")
    plt.semilogx(
        np.squeeze(d), np.squeeze(loss_ter),
        label="free space + clutter loss",
    )

    plt.title("Free space with additional median clutter loss ($f=27GHz$)")
    plt.xlabel("distance [m]")
    plt.ylabel("path loss [dB]")
    plt.xlim((0, d[-1]))
    plt.ylim((80, 240))
    plt.legend(loc="upper left")
    plt.tight_layout()
    plt.grid()

    plt.show()
