# -*- coding: utf-8 -*-
"""
Created on Tue Jul  4 11:57:41 2017

@author: LeticiaValle_Mac
"""

from sharc.propagation.propagation import Propagation

import numpy as np


class PropagationABG(Propagation):
    """
    Implements the ABG loss model according to the article "Propagation Path
    Loss Models for 5G Urban Microand Macro-Cellular Scenarios"
    """

    def __init__(self, random_number_gen: np.random.RandomState):
        super().__init__(random_number_gen)
        self.alpha = 3.4
        self.beta = 19.2
        self.gamma = 2.3
        self.building_loss = 20
        self.shadowing_sigma_dB = 6.5


    def get_loss(self, *args, **kwargs) -> np.array:
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
        f = kwargs["frequency"]
        indoor_stations = kwargs["indoor_stations"]

        if "distance_3D" in kwargs:
            d = kwargs["distance_3D"]
        else:
            d = kwargs["distance_2D"]

        if "alpha" in kwargs:
            self.alpha = kwargs["alpha"]

        if "beta" in kwargs:
            self.beta = kwargs["beta"]

        if "gamma" in kwargs:
            self.gamma = kwargs["gamma"]

        if "shadowing" in kwargs:
            std = kwargs["shadowing"]
        else:
            std = False

        if std:
            shadowing = self.random_number_gen.normal(0, self.shadowing_sigma_dB, d.shape)
        else:
            shadowing = 0

        building_loss = self.building_loss*indoor_stations

        loss = 10*self.alpha*np.log10(d) + self.beta + 10*self.gamma*np.log10(f*1e-3) + \
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
    distance_2D = np.linspace(1, 1000, num=1000)[:,np.newaxis]
    freq = 26000*np.ones(distance_2D.shape)
    h_bs = 25*np.ones(len(distance_2D[:,0]))
    h_ue = 1.5*np.ones(len(distance_2D[0,:]))
    h_e = np.zeros(distance_2D.shape)
    distance_3D = np.sqrt(distance_2D**2 + (h_bs[:,np.newaxis] - h_ue)**2)

    uma = PropagationUMa()
    umi = PropagationUMi()
    abg = PropagationABG()
    freespace = PropagationFreeSpace()

    uma_los = uma.get_loss_los(distance_2D, distance_3D, freq, h_bs, h_ue, h_e, shadowing_std)
    uma_nlos = uma.get_loss_nlos(distance_2D, distance_3D, freq, h_bs, h_ue, h_e, shadowing_std)
    umi_los = umi.get_loss_los(distance_2D, distance_3D, freq, h_bs, h_ue, h_e, shadowing_std)
    umi_nlos = umi.get_loss_nlos(distance_2D, distance_3D, freq, h_bs, h_ue, h_e, shadowing_std)
    fs = freespace.get_loss(distance_2D=distance_2D, frequency=freq)
    abg_los = abg.get_loss(distance_2D=distance_2D, frequency=freq, indoor_stations=0)

    fig = plt.figure(figsize=(8,6), facecolor='w', edgecolor='k')
    ax = fig.gca()
    #ax.set_prop_cycle( cycler('color', ['r', 'g', 'b', 'y']) )

    ax.semilogx(distance_2D, uma_los, "-r", label="UMa LOS")
    ax.semilogx(distance_2D, uma_nlos, "--r", label="UMa NLOS")
    ax.semilogx(distance_2D, umi_los, "-b", label="UMi LOS")
    ax.semilogx(distance_2D, umi_nlos, "--b", label="UMi NLOS")
    ax.semilogx(distance_2D, abg_los, "-g", label="ABG")
    ax.semilogx(distance_2D, fs, "-k", label="free space")

    plt.title("Path loss models")
    plt.xlabel("distance [m]")
    plt.ylabel("path loss [dB]")
    plt.xlim((0, distance_2D[-1,0]))
    #plt.ylim((0, 1.1))
    plt.legend(loc="upper left")
    plt.tight_layout()
    plt.grid()

    plt.show()
