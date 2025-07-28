# -*- coding: utf-8 -*-
"""
Created on Mon Jul  3 10:29:47 2017

@author: LeticiaValle_Mac
"""
import numpy as np
from multipledispatch import dispatch

from sharc.propagation.propagation import Propagation
from sharc.station_manager import StationManager
from sharc.parameters.parameters import Parameters


class PropagationUMi(Propagation):
    """
    Implements the Urban Micro path loss model (Street Canyon) with LOS
    probability according to 3GPP TR 38.900 v14.2.0.
    TODO: calculate the effective environment height for the generic case
    """

    def __init__(
        self,
        random_number_gen: np.random.RandomState,
        los_adjustment_factor: float,
    ):
        super().__init__(random_number_gen)
        self.los_adjustment_factor = los_adjustment_factor

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
        """Wrapper function for the PropagationUMi get_loss method
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

        if wrap_around_enabled and (station_a.is_imt_station() and station_b.is_imt_station()):
            distance_2d, distance_3d, _, _ = \
                station_a.get_dist_angles_wrap_around(station_b)
        else:
            distance_2d = station_a.get_distance_to(station_b)
            distance_3d = station_a.get_3d_distance_to(station_b)

        loss = self.get_loss(
            distance_3d,
            distance_2d,
            frequency * np.ones(distance_2d.shape),
            station_b.height,
            station_a.height,
            params.imt.shadowing,
        )

        return loss

    # pylint: disable=function-redefined
    # pylint: disable=arguments-renamed
    @dispatch(np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, bool)
    def get_loss(
        self,
        distance_3D: np.array,
        distance_2D: np.array,
        frequency: np.array,
        bs_height: np.array,
        ue_height: np.array,
        shadowing_flag: bool,
    ) -> np.array:
        """
        Calculates path loss for LOS and NLOS cases with respective shadowing
        (if shadowing is to be added)

        Parameters
        ----------
            distance_3D (np.array) : 3D distances between base stations and user equipment
            distance_2D (np.array) : 2D distances between base stations and user equipment
            frequency (np.array) : center frequencies [MHz]
            bs_height (np.array) : base station antenna heights
            ue_height (np.array) : user equipment antenna heights
            shadowing (bool) : if shadowing should be added or not

        Returns
        -------
            array with path loss values with dimensions of distance_2D
        """
        if shadowing_flag:
            shadowing_los = 4
            shadowing_nlos = 7.82    # option 1 for UMi NLOS
            # shadowing_nlos = 8.2    # option 2 for UMi NLOS
        else:
            shadowing_los = 0
            shadowing_nlos = 0

        # effective height
        h_e = np.ones(distance_2D.shape)

        los_probability = self.get_los_probability(
            distance_2D, self.los_adjustment_factor,
        )
        los_condition = self.get_los_condition(los_probability)

        i_los = np.where(los_condition == True)[:2]
        i_nlos = np.where(los_condition == False)[:2]

        loss = np.empty(distance_2D.shape)

        if len(i_los[0]):
            loss_los = self.get_loss_los(
                distance_2D, distance_3D, frequency, bs_height, ue_height, h_e, shadowing_los,
            )
            loss[i_los] = loss_los[i_los]

        if len(i_nlos[0]):
            loss_nlos = self.get_loss_nlos(
                distance_2D, distance_3D, frequency, bs_height, ue_height, h_e, shadowing_nlos,
            )
            loss[i_nlos] = loss_nlos[i_nlos]

        return loss

    def get_loss_los(
        self, distance_2D: np.array,
        distance_3D: np.array,
        frequency: np.array,
        h_bs: np.array,
        h_ue: np.array,
        h_e: np.array,
        shadowing_std=4,
    ):
        """
        Calculates path loss for the LOS (line-of-sight) case.

        Parameters
        ----------
            distance_2D : array of 2D distances between base stations and user
                          equipments [m]
            distance_3D : array of 3D distances between base stations and user
                          equipments [m]
            frequency : center frequency [MHz]
            h_bs : array of base stations antenna heights [m]
            h_ue : array of user equipments antenna heights [m]
        """
        breakpoint_distance = self.get_breakpoint_distance(
            frequency, h_bs, h_ue, h_e,
        )

        # get index where distance if less than breakpoint
        idl = np.where(distance_2D < breakpoint_distance)

        # get index where distance if greater than breakpoint
        idg = np.where(distance_2D >= breakpoint_distance)

        loss = np.empty(distance_2D.shape)

        if len(idl[0]):
            loss[idl] = 21 * np.log10(distance_3D[idl]) + \
                20 * np.log10(frequency[idl]) - 27.55

        if len(idg[0]):
            fitting_term = -9.5 * \
                np.log10(
                    breakpoint_distance**2 +
                    (h_bs - h_ue[:, np.newaxis])**2,
                )
            loss[idg] = 40 * np.log10(distance_3D[idg]) + 20 * np.log10(frequency[idg]) - 27.55 \
                + fitting_term[idg]

        if shadowing_std:
            shadowing = self.random_number_gen.normal(
                0, shadowing_std, distance_2D.shape,
            )
        else:
            shadowing = 0

        return loss + shadowing

    def get_loss_nlos(
        self, distance_2D: np.array, distance_3D: np.array,
        frequency: np.array,
        h_bs: np.array, h_ue: np.array, h_e: np.array,
        shadowing_std=7.82,
    ):
        """
        Calculates path loss for the NLOS (non line-of-sight) case.

        Parameters
        ----------
            distance_2D : array of 2D distances between base stations and user
                          equipments [m]
            distance_3D : array of 3D distances between base stations and user
                          equipments [m]
            frequency : center frequency [MHz]
            h_bs : array of base stations antenna heights [m]
            h_ue : array of user equipments antenna heights [m]
        """
        # option 1 for UMi NLOS
        loss_nlos = -37.55 + 35.3 * np.log10(distance_3D) + 21.3 * np.log10(frequency) \
            - 0.3 * (h_ue[:, np.newaxis] - 1.5)

        loss_los = self.get_loss_los(
            distance_2D, distance_3D, frequency, h_bs, h_ue, h_e, 0,
        )
        loss_nlos = np.maximum(loss_los, loss_nlos)

        # option 2 for UMi NLOS
        # loss_nlos = 31.9*np.log10(distance_3D) + 20*np.log10(frequency*1e-3) + 32.4

        if shadowing_std:
            shadowing = self.random_number_gen.normal(
                0, shadowing_std, distance_3D.shape,
            )
        else:
            shadowing = 0

        return loss_nlos + shadowing

    def get_breakpoint_distance(self, frequency: float, h_bs: np.array, h_ue: np.array, h_e: np.array) -> float:
        """
        Calculates the breakpoint distance for the LOS (line-of-sight) case.

        Parameters
        ----------
            frequency : centre frequency [MHz]
            h_bs : array of actual base station antenna height [m]
            h_ue : array of actual user equipment antenna height [m]
            h_e : array of effective environment height [m]

        Returns
        -------
            array of breakpoint distances [m]
        """
        #  calculate the effective antenna heights
        h_bs_eff = h_bs - h_e
        h_ue_eff = h_ue[:, np.newaxis] - h_e

        # calculate the breakpoint distance
        breakpoint_distance = 4 * h_bs_eff * \
            h_ue_eff * (frequency * 1e6) / (3e8)
        return breakpoint_distance

    def get_los_condition(self, p_los: np.array) -> np.array:
        """
        Evaluates if user equipments are LOS (True) of NLOS (False).

        Parameters
        ----------
            p_los : array with LOS probabilities for each user equipment.

        Returns
        -------
            An array with True or False if user equipments are in LOS of NLOS
            condition, respectively.
        """
        los_condition = self.random_number_gen.random_sample(
            p_los.shape,
        ) < p_los
        return los_condition

    def get_los_probability(
        self,
        distance_2D: np.array,
        los_adjustment_factor: float,
    ) -> np.array:
        """
        Returns the line-of-sight (LOS) probability

        Parameters
        ----------
            distance_2D : Two-dimensional array with 2D distance values from
                          base station to user terminal [m]
            los_adjustment_factor : adjustment factor to increase/decrease the
                          LOS probability. Original value is 18 as per 3GPP

        Returns
        -------
            LOS probability as a numpy array with same length as distance
        """

        p_los = np.ones(distance_2D.shape)
        idl = np.where(distance_2D > los_adjustment_factor)
        p_los[idl] = (
            los_adjustment_factor / distance_2D[idl] +
            np.exp(-distance_2D[idl] / 36) * (1 - los_adjustment_factor / distance_2D[idl])
        )

        return p_los


if __name__ == '__main__':

    import matplotlib.pyplot as plt
    from cycler import cycler

    ###########################################################################
    # Print LOS probability
    # h_ue = np.array([1.5, 17, 23])
    distance_2D = np.linspace(1, 1000, num=1000)[:, np.newaxis]
    umi = PropagationUMi(np.random.RandomState(101), 18)

    los_probability = np.empty(distance_2D.shape)
    name = list()

    los_adjustment_factor = 18
    los_probability_18 = umi.get_los_probability(
        distance_2D, los_adjustment_factor,
    )

    los_adjustment_factor = 29
    los_probability_29 = umi.get_los_probability(
        distance_2D, los_adjustment_factor,
    )

    fig = plt.figure(figsize=(6, 5), facecolor='w', edgecolor='k')
    ax = fig.gca()
    ax.set_prop_cycle(cycler('color', ['r', 'g', 'b', 'y']))

    ax.loglog(distance_2D, 100 * los_probability_18, label=r"$\alpha = 18$")
    ax.loglog(distance_2D, 100 * los_probability_29, label=r"$\alpha = 29$")

    plt.title("UMi - LOS probability")
    plt.xlabel("distance between BS and UE [m]")
    plt.ylabel("probability [%]")
    ax.legend(loc="lower left")
    ax.grid(True, which="both", color="grey", linestyle='-', linewidth=0.2)
    plt.xlim((10, 300))
    plt.ylim((10, 102))
    plt.tight_layout()
    plt.grid()

    ###########################################################################
    # Print path loss for UMi-LOS, UMi-NLOS and Free Space
    from propagation_free_space import PropagationFreeSpace
    shadowing_std = 0
    #  1 ue x 1000 bs
    num_ue = 1
    num_bs = 1000
    distance_2D = np.repeat(np.linspace(1, num_bs, num=num_bs)[np.newaxis, :], num_ue, axis=0)
    freq = 24350 * np.ones(distance_2D.shape)
    h_bs = 6 * np.ones(num_bs)
    h_ue = 1.5 * np.ones(num_ue)
    h_e = np.ones(distance_2D.shape)
    distance_3D = np.sqrt(distance_2D**2 + (h_bs - h_ue[np.newaxis, :])**2)

    loss_los = umi.get_loss_los(
        distance_2D, distance_3D, freq, h_bs, h_ue, h_e, shadowing_std,
    )
    loss_nlos = umi.get_loss_nlos(
        distance_2D, distance_3D, freq, h_bs, h_ue, h_e, shadowing_std,
    )
    loss_fs = PropagationFreeSpace(
        np.random.RandomState(101,),
    ).get_free_space_loss(distance=distance_2D, frequency=freq)

    fig = plt.figure(figsize=(8, 6), facecolor='w', edgecolor='k')
    ax = fig.gca()
    ax.set_prop_cycle(cycler('color', ['r', 'g', 'b', 'y']))

    ax.semilogx(distance_2D[0, :], loss_los[0, :], label="UMi LOS")
    ax.semilogx(distance_2D[0, :], loss_nlos[0, :], label="UMi NLOS")
    ax.semilogx(distance_2D[0, :], loss_fs[0, :], label="free space")

    plt.title("UMi - path loss")
    plt.xlabel("distance [m]")
    plt.ylabel("path loss [dB]")
    plt.xlim((0, distance_2D[0, -1]))
    plt.ylim((60, 200))
    plt.legend(loc="upper left")
    plt.tight_layout()
    plt.grid()

    plt.show()
