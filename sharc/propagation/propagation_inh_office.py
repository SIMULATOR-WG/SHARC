# -*- coding: utf-8 -*-
"""
Created on Wed Dec  6 19:23:21 2017

@author: edgar
"""

from sharc.propagation.propagation import Propagation

import numpy as np
import matplotlib.pyplot as plt
from cycler import cycler

class PropagationInhOffice(Propagation):
    """
    Implements the Indoor Hotspot - Office path loss model with LOS probability
    according to 3GPP TR 38.900 v14.2.0.
    """


    def get_loss(self, *args, **kwargs) -> np.array:
        """
        Calculates path loss for LOS and NLOS cases with respective shadowing
        (if shadowing has to be added)

        Parameters
        ----------
            distance_3D (np.array) : 3D distances between stations
            distance_2D (np.array) : 2D distances between stations
            frequency (np.array) : center frequencie [MHz]
            indoor (np.array) : indicates whether UE is indoor
            shadowing (bool) : if shadowing should be added or not

        Returns
        -------
            array with path loss values with dimensions of distance_2D

        """
        d_3D = kwargs["distance_3D"]
        d_2D = kwargs["distance_2D"]
        f = kwargs["frequency"]
        indoor = kwargs["indoor"]
        std = kwargs["shadowing"]

        if std:
            shadowing_los = 3
            shadowing_nlos = 8.03
        else:
            shadowing_los = 0
            shadowing_nlos = 0

        los_probability = self.get_los_probability(d_2D)
        los_condition = self.get_los_condition(los_probability, indoor)

        i_los = np.where(los_condition)[:2]
        i_nlos = np.where(~ los_condition)[:2]

        loss = np.empty(d_2D.shape)

        if len(i_los[0]):
            loss[i_los] = self.get_loss_los(d_3D[i_los], f[i_los], shadowing_los)

        if len(i_nlos[0]):
            loss[i_nlos] = self.get_loss_nlos(d_3D[i_nlos], f[i_nlos], shadowing_nlos)

        return loss


    def get_loss_los(self,
                     distance_3D: np.array,
                     frequency: np.array,
                     shadowing_std: float):
        """
        Calculates path loss for the LOS (line-of-sight) case.

        Parameters
        ----------
            distance_3D : array of 3D distances between BS and UE [m]
            frequency : center frequency [MHz]
            shadowing_std : standard deviation of shadowing [dB]

        Returns
        -------
            array with path loss values with dimensions of distance_3D
        """

        loss = 32.4 + 17.3*np.log10(distance_3D) + 20*np.log10(frequency/1e3)

        if shadowing_std:
            shadowing = self.random_number_gen.normal(0, shadowing_std, distance_3D.shape)
        else:
            shadowing = 0

        return loss + shadowing


    def get_loss_nlos(self,
                      distance_3D: np.array,
                      frequency: np.array,
                      shadowing_std: float):
        """
        Calculates path loss for the NLOS (non line-of-sight) case.

        Parameters
        ----------
            distance_3D : array of 3D distances between BS and UE [m]
            frequency : center frequency [MHz]
            shadowing_std : standard deviation of shadowing [dB]

        Returns
        -------
            array with path loss values with dimensions of distance_3D
        """

        loss_nlos = 17.3 + 38.3*np.log10(distance_3D) + 24.9*np.log10(frequency/1e3)
        loss_los = self.get_loss_los(distance_3D, frequency, 0)

        loss = np.maximum(loss_los, loss_nlos)

        if shadowing_std:
            shadowing = self.random_number_gen.normal(0, shadowing_std, distance_3D.shape)
        else:
            shadowing = 0

        return loss + shadowing


    def get_los_condition(self, p_los: np.array, indoor: np.array) -> np.array:
        """
        Evaluates if user equipments are LOS (True) of NLOS (False). If UE is
        outdoor, condition is automatically set to NLOS.

        Parameters
        ----------
            p_los : array with LOS probabilities for each user equipment.
            indoor : indicates whether UE is indoor

        Returns
        -------
            An array with True or False if user equipments are in LOS of NLOS
            condition, respectively.
        """
        los_condition = self.random_number_gen.random_sample(p_los.shape) < p_los
        los_condition = los_condition & indoor
        return los_condition


    def get_los_probability(self, distance_2D: np.array) -> np.array:
        """
        Returns the line-of-sight (LOS) probability

        Parameters
        ----------
            distance_2D : Two-dimensional array with 2D distance values from
                          base station to user terminal [m]

        Returns
        -------
            LOS probability as a numpy array with same length as distance
        """

        p_los = np.ones(distance_2D.shape)
        id1 = np.where((distance_2D > 1.2) & (distance_2D < 6.5))
        p_los[id1] = np.exp(-(distance_2D[id1] - 1.2)/4.7)
        id2 = np.where(distance_2D >= 6.5)
        p_los[id2] = np.exp(-(distance_2D[id2] - 6.5)/32.6)*0.32

        return p_los


if __name__ == '__main__':

    ###########################################################################
    # Print LOS probability
    distance_2D = np.linspace(0.1, 150, num=1000)[:,np.newaxis]
    inh = PropagationInhOffice(np.random.RandomState())

    los_probability = inh.get_los_probability(distance_2D)

    plt.figure(figsize=(8,6), facecolor='w', edgecolor='k')
    plt.loglog(distance_2D, los_probability)

    plt.title("InH Office - LOS probability")
    plt.xlabel("distance [m]")
    plt.ylabel("probability")
    plt.xlim((0, distance_2D[-1,0]))
    plt.ylim((0, 1.1))
    plt.tight_layout()
    plt.grid()
    plt.show()

    ###########################################################################
    # Print path loss for InH-LOS, InH-NLOS and Free Space
    from sharc.propagation.propagation_free_space import PropagationFreeSpace
    shadowing = 0
    distance_2D = np.linspace(1, 150, num=1000)[:,np.newaxis]
    frequency = 27000*np.ones(distance_2D.shape)
    h_bs = 3*np.ones(len(distance_2D[:,0]))
    h_ue = 1.5*np.ones(len(distance_2D[0,:]))
    distance_3D = np.sqrt(distance_2D**2 + (h_bs[:,np.newaxis] - h_ue)**2)
    indoor = np.ones(len(distance_2D[0,:]), dtype = bool)

    loss_fs = PropagationFreeSpace(np.random.RandomState()).get_loss(distance_3D=distance_3D, frequency=frequency)

    loss_inh = inh.get_loss(distance_3D = distance_3D,
                            distance_2D = distance_2D,
                            frequency = frequency,
                            indoor = indoor,
                            shadowing = shadowing)

    fig = plt.figure(figsize=(8,6), facecolor='w', edgecolor='k')
    ax = fig.gca()
    ax.set_prop_cycle( cycler('color', ['r', 'b', 'g', 'y']) )

    ax.scatter(distance_2D, loss_inh, label = "InH-Office")
    ax.plot(distance_2D, loss_fs, label = "free space")
    ax.set_xscale('log')

    plt.title("InH-Office - path loss")
    plt.xlabel("distance [m]")
    plt.ylabel("path loss [dB]")
    plt.xlim((0, distance_2D[-1,0]))

    plt.legend(loc="upper left")
    plt.tight_layout()
    plt.grid()

    plt.show()
