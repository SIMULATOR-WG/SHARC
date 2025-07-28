# -*- coding: utf-8 -*-
"""
Created on Tue Jul  4 11:57:41 2017

@author: LeticiaValle_Mac
"""

import numpy as np
from multipledispatch import dispatch

from sharc.propagation.propagation import Propagation
from sharc.station_manager import StationManager
from sharc.parameters.parameters import Parameters


class PropagationABG(Propagation):
    """
    Implements the ABG loss model according to the article "Propagation Path
    Loss Models for 5G Urban Microand Macro-Cellular Scenarios"
    """

    def __init__(
        self, random_number_gen: np.random.RandomState,
        alpha=3.4, beta=19.2, gamma=2.3, building_loss=20, shadowing_sigma_dB=6.5,
    ):
        super().__init__(random_number_gen)
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma
        self.building_loss = 20
        self.shadowing_sigma_dB = 6.5

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
        """Wrapper function for the get_loss method
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
        wrap_around_enabled = False
        if params.imt.topology.type == "MACROCELL":
            wrap_around_enabled = params.imt.topology.macrocell.wrap_around \
                                    and params.imt.topology.macrocell.num_clusters == 1
        if params.imt.topology.type == "HOTSPOT":
            wrap_around_enabled = params.imt.topology.hotspot.wrap_around \
                                    and params.imt.topology.hotspot.num_clusters == 1

        if wrap_around_enabled:
            _, distances_3d, _, _ = \
                station_a.get_dist_angles_wrap_around(station_b)
        else:
            distances_3d = station_a.get_3d_distance_to(station_b)

        indoor_stations = station_a.indoor

        loss = \
            self.get_loss(
                distances_3d,
                frequency * np.ones(distances_3d.shape),
                indoor_stations,
                params.imt.shadowing,
            )

        return loss

    @dispatch(np.ndarray, np.ndarray, np.ndarray, bool)
    def get_loss(self, distance: np.array, frequency: np.array, indoor_stations: np.array, shadowing: bool) -> np.array:
        """
        Calculates path loss for LOS and NLOS cases with respective shadowing
        (if shadowing is to be added)

        Parameters
        ----------
            distance_2D (np.array) : distances between stations [m]
            frequency (np.array) : center frequencie [MHz]
            indoor_stations (np.array) : array indicating stations that are indoor
            alpha (float): captures how the PL increases as the distance increases
            beta (float): floating offset value in dB
            gamma(float): captures the PL variation over the frequency
            shadowing (bool) : standard deviation value

        Returns
        -------
            array with path loss values with dimensions of distance_2D

        """
        if shadowing:
            shadowing = self.random_number_gen.normal(
                0, self.shadowing_sigma_dB, distance.shape,
            )
        else:
            shadowing = 0

        building_loss = self.building_loss * np.tile(indoor_stations, (distance.shape[1], 1)).transpose()

        loss = 10 * self.alpha * np.log10(distance) + self.beta + 10 * self.gamma * np.log10(frequency * 1e-3) + \
            shadowing + building_loss

        return loss


if __name__ == '__main__':

    ###########################################################################
    # Print path loss for ABG and Free Space models
    from sharc.propagation.propagation_free_space import PropagationFreeSpace
    from sharc.propagation.propagation_uma import PropagationUMa
    from sharc.propagation.propagation_umi import PropagationUMi

    import matplotlib.pyplot as plt

    shadowing_std = 0
    distance_2D = np.linspace(1, 1000, num=1000)[:, np.newaxis]
    freq = 26000 * np.ones(distance_2D.shape)
    h_bs = 25 * np.ones(len(distance_2D[:, 0]))
    h_ue = 1.5 * np.ones(len(distance_2D[0, :]))
    h_e = np.zeros(distance_2D.shape)
    distance_3D = np.sqrt(distance_2D**2 + (h_bs[:, np.newaxis] - h_ue)**2)

    random_number_gen = np.random.RandomState(101)

    uma = PropagationUMa(random_number_gen)
    umi = PropagationUMi(random_number_gen, 18)
    abg = PropagationABG(random_number_gen)
    freespace = PropagationFreeSpace(random_number_gen)

    uma_los = uma.get_loss_los(
        distance_2D, distance_3D, freq, h_bs, h_ue, h_e, shadowing_std,
    )
    uma_nlos = uma.get_loss_nlos(
        distance_2D, distance_3D, freq, h_bs, h_ue, h_e, shadowing_std,
    )
    umi_los = umi.get_loss_los(
        distance_2D, distance_3D, freq, h_bs, h_ue, h_e, shadowing_std,
    )
    umi_nlos = umi.get_loss_nlos(
        distance_2D, distance_3D, freq, h_bs, h_ue, h_e, shadowing_std,
    )
    fs = freespace.get_loss(distance_2D, freq)
    abg_los = abg.get_loss(distance_2D, freq, np.zeros(shape=distance_2D.shape, dtype=bool,), False,)

    fig = plt.figure(figsize=(8, 6), facecolor='w', edgecolor='k')
    ax = fig.gca()
    # ax.set_prop_cycle( cycler('color', ['r', 'g', 'b', 'y']) )

    ax.semilogx(distance_2D, uma_los, "-r", label="UMa LOS")
    ax.semilogx(distance_2D, uma_nlos, "--r", label="UMa NLOS")
    ax.semilogx(distance_2D, umi_los, "-b", label="UMi LOS")
    ax.semilogx(distance_2D, umi_nlos, "--b", label="UMi NLOS")
    ax.semilogx(distance_2D, abg_los, "-g", label="ABG")
    ax.semilogx(distance_2D, fs, "-k", label="free space")

    plt.title("Path loss models")
    plt.xlabel("distance [m]")
    plt.ylabel("path loss [dB]")
    plt.xlim((0, distance_2D[-1, 0]))
    # plt.ylim((0, 1.1))
    plt.legend(loc="upper left")
    plt.tight_layout()
    plt.grid()

    plt.show()
