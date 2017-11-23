# -*- coding: utf-8 -*-
"""
Created on Wed Aug 16 13:42:19 2017

@author: edgar
"""

from sharc.propagation.propagation import Propagation
from sharc.propagation.propagation_free_space import PropagationFreeSpace
from sharc.propagation.propagation_clutter_loss import PropagationClutterLoss
from sharc.support.enumerations import StationType

import numpy as np

class PropagationTerSimple(Propagation):
    """
    Implements the simplified terrestrial propagation model, which is the
    basic free space and additional clutter losses
    """

    def __init__(self):
        super().__init__()
        self.clutter = PropagationClutterLoss()
        self.free_space = PropagationFreeSpace()
        self.building_loss = 20


    def get_loss(self, *args, **kwargs) -> np.array:
        if "distance_2D" in kwargs:
            d = kwargs["distance_2D"]
        else:
            d = kwargs["distance_3D"]

        f = kwargs["frequency"]
        p = kwargs.pop("loc_percentage", "RANDOM")
        indoor_stations = kwargs["indoor_stations"]
        number_of_sectors = kwargs.pop("number_of_sectors",1)

        free_space_loss = self.free_space.get_loss(distance_2D=d,
                                                   frequency=f)

        clutter_loss = self.clutter.get_loss(frequency=f,
                                             distance=d,
                                             loc_percentage=p,
                                             station_type=StationType.FSS_ES)

        building_loss = self.building_loss*indoor_stations

        loss = free_space_loss + building_loss + clutter_loss
        loss = np.repeat(loss, number_of_sectors, 1)

        return loss


if __name__ == '__main__':

    import matplotlib.pyplot as plt

    ###########################################################################
    # Print path loss for TerrestrialSimple and Free Space

    d = np.linspace(10, 10000, num=10000)
    freq = 27000*np.ones(d.shape)
    indoor_stations = np.zeros(d.shape, dtype = bool)
    loc_percentage = 0.5

    free_space = PropagationFreeSpace()
    ter_simple = PropagationTerSimple()

    loss_ter = ter_simple.get_loss(distance_2D = d,
                                  frequency = freq,
                                  loc_percentage = loc_percentage,
                                  indoor_stations = indoor_stations)

    loss_fs = free_space.get_loss(distance_2D = d,
                                  frequency = freq)

    fig = plt.figure(figsize=(8,6), facecolor='w', edgecolor='k')

    plt.semilogx(np.squeeze(d), np.squeeze(loss_fs), label = "free space")
    plt.semilogx(np.squeeze(d), np.squeeze(loss_ter), label = "free space + clutter loss")

    plt.title("Free space with additional median clutter loss ($f=27GHz$)")
    plt.xlabel("distance [m]")
    plt.ylabel("path loss [dB]")
    plt.xlim((0, d[-1]))
    plt.ylim((80, 240))
    plt.legend(loc="upper left")
    plt.tight_layout()
    plt.grid()

    plt.show()
