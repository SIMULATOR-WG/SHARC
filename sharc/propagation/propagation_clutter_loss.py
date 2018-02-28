# -*- coding: utf-8 -*-
"""
Created on Mon May 15 12:51:48 2017

@author: LeticiaValle_Mac
"""

from sharc.propagation.propagation import Propagation
from sharc.support.enumerations import StationType

import numpy as np
import scipy
import math
from scipy import special


class PropagationClutterLoss(Propagation):


    def get_loss(self, *args, **kwargs) -> np.array:
        """
        Calculates clutter loss.

        Parameters
        ----------
            distance (np.array) : distances between stations [m]
            frequency (np.array) : center frequency [MHz]
            elevation (np.array) : elevation angles [deg]
            loc_percentage (np.array) : Percentage locations range [0, 1[
                                        "RANDOM" for random percentage (Default = RANDOM)
            station_type (StationType) : if type is IMT or FSS_ES, assume terrestrial
                terminal within the clutter (ref § 3.2); otherwise, assume that
                one terminal is within the clutter and the other is a satellite,
                aeroplane or other platform above the surface of the Earth.

        Returns
        -------
            array with clutter loss values with dimensions of distance

        """
        f = kwargs["frequency"]
        loc_per = kwargs.pop("loc_percentage","RANDOM")
        type = kwargs["station_type"]

        d = kwargs["distance"]

        if f.size == 1:
            f = f * np.ones(d.shape)

        if isinstance(loc_per, str) and loc_per.upper() == "RANDOM":
            p = self.random_number_gen.random_sample(d.shape)
        else:
            p = loc_per*np.ones(d.shape)

        if type is StationType.IMT_BS or type is StationType.IMT_UE or type is StationType.FSS_ES:
            loss = self.get_terrestrial_clutter_loss(f, d, p)
        else:
            theta = kwargs["elevation"]
            loss = self.get_spacial_clutter_loss(f, theta, p)
        return loss

    def get_spacial_clutter_loss(self, frequency : float,
                             elevation_angle : float,
                             loc_percentage):
        """
        This method models the calculation of the statistical distribution of
        clutter loss where one end of the interference path is within man-made
        clutter, and the other is a satellite, aeroplane, or other platform
        above the surface of the Earth. This model is applicable to urban and
        suburban environments.

        Parameters
        ----
            frequency : center frequency [MHz]
            elevation_angle : elevation angle [degrees]
            loc_percentage : percentage of locations [0,1[

        Returns
        -------
            loss : The clutter loss not exceeded for p% of locations for the
                terrestrial to terrestrial path
        """
        k1 = 93*(frequency*1e-3)**0.175
        A1 = 0.05

        y = np.sin(A1*(1 - (elevation_angle/90)) + math.pi*(elevation_angle/180))
        y1 = np.cos(A1*(1 - (elevation_angle/90)) + math.pi*(elevation_angle/180))

        cot = (y1/y)
        Q = np.sqrt(2)*scipy.special.erfcinv(2*loc_percentage)
        loss = (-k1*(np.log(1 - loc_percentage))*cot)**(0.5*(90 - elevation_angle)/90) - 1 - 0.6*Q

        return loss


    def get_terrestrial_clutter_loss(self,
                                     frequency: float,
                                     distance: float,
                                     loc_percentage: float,
                                     apply_both_ends = True):
        """
        This method gives models the statistical distribution of clutter loss.
        The model can be applied for urban and suburban clutter loss modelling.

        Parameters
        ----
            frequency : center frequency [MHz]
            distance : distance [m]
            loc_percentage : percentage of locations [0,1]
            apply_both_ends : if correction will be applied at both ends of the path

        Returns
        -------
            loss : The clutter loss not exceeded for p% of locations for the
                terrestrial to terrestrial path
        """

        d = distance.reshape((-1, 1))
        f = frequency.reshape((-1, 1))
        p = loc_percentage.reshape((-1, 1))

        loss = np.zeros(d.shape)

        # minimum path length for the correction to be applied at only one end of the path
        id_1 = np.where(d >= 250)[0]

        if len(id_1):
            Lt = 23.5 + 9.6 * np.log10(f[id_1] * 1e-3)
            Ls = 32.98 + 23.9 * np.log10(d[id_1] * 1e-3) + 3 * np.log10(f[id_1] * 1e-3)
            Q = np.sqrt(2) * scipy.special.erfcinv(2 * (p[id_1]))
            loss[id_1] = -5 * np.log10(10 ** (-0.2 * Lt) + 10 ** (-0.2 * Ls)) - 6 * Q

        # minimum path length for the correction to be applied at only one end of the path
        id_2 = np.where(d >= 1000)[0]

        if apply_both_ends and len(id_2):
            Lt = 23.5 + 9.6 * np.log10(f[id_2] * 1e-3)
            Ls = 32.98 + 23.9 * np.log10(d[id_2] * 1e-3) + 3 * np.log10(f[id_2] * 1e-3)
            Q = np.sqrt(2) * scipy.special.erfcinv(2 * (p[id_2]))
            loss[id_2] = loss[id_2] + (-5 * np.log10(10 ** (-0.2 * Lt) + 10 ** (-0.2 * Ls)) - 6 * Q)

        loss = loss.reshape(distance.shape)

        return loss


if __name__ == '__main__':
    import matplotlib.pyplot as plt

    elevation_angle = np.array([90, 80, 70, 60, 50, 40, 30, 20, 15, 10, 5, 0])
    loc_percentage = np.linspace(0, 1, 1001)
    frequency = 27250 * np.ones(elevation_angle.shape)

    cl = PropagationClutterLoss()
    clutter_loss = np.empty([len(elevation_angle), len(loc_percentage)])

    for i in range(len(loc_percentage)):
        clutter_loss[:, i] = cl.get_spacial_clutter_loss(frequency,
                                                         elevation_angle,
                                                         loc_percentage[i])

    fig = plt.figure(figsize=(8, 6), facecolor='w', edgecolor='k')
    ax = fig.gca()

    for j in range(len(elevation_angle)):
        ax.plot(clutter_loss[j, :], 100 * loc_percentage, label="%i deg" % elevation_angle[j], linewidth=1)

    plt.title("Cumulative distribution of clutter loss not exceeded for 27 GHz")
    plt.xlabel("clutter loss [dB]")

    plt.ylabel("percent of locations [%]")
    plt.xlim((-10, 70))
    plt.ylim((0, 100))
    plt.legend(loc="lower right")
    plt.tight_layout()
    plt.grid()

    distance = np.linspace(250, 100000, 100000)
    frequency = np.array([2, 3, 6, 16, 40, 67]) * 1e3

    loc_percentage = 0.5 * np.ones(distance.shape)
    apply_both_ends = False

    clutter_loss_ter = np.empty([len(frequency), len(distance)])

    for i in range(len(frequency)):
        clutter_loss_ter[i, :] = cl.get_terrestrial_clutter_loss(frequency[i] * np.ones(distance.shape),
                                                                 distance,
                                                                 loc_percentage,
                                                                 apply_both_ends)

    fig = plt.figure(figsize=(8, 6), facecolor='w', edgecolor='k')
    ax = fig.gca()
    # ax.set_prop_cycle( cycler('color', ['k', 'r', 'b', 'g']) )

    for j in range(len(frequency)):
        freq = frequency[j] * 1e-3
        ax.semilogx(distance * 1e-3, clutter_loss_ter[j, :], label="%i GHz" % freq, linewidth=1)

    plt.title("Median clutter loss for terrestrial paths")
    plt.xlabel("Distance [km]")
