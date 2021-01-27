# -*- coding: utf-8 -*-
"""
Created on Thu Feb 16 12:04:27 2017

@author: edgar
"""

from sharc.propagation.propagation import Propagation

import numpy as np

class PropagationFreeSpace(Propagation):
    """
    Implements the Free Space propagation model.
    Frequency in MHz and distance in meters
    """

    def get_loss(self, *args, **kwargs) -> np.array:
        if "distance_2D" in kwargs:
            d = kwargs["distance_2D"]
        else:
            d = kwargs["distance_3D"]

        f = kwargs["frequency"]
        number_of_sectors = kwargs.pop("number_of_sectors",1)

        loss = 20*np.log10(d) + 20*np.log10(f) - 27.55

        if number_of_sectors > 1:
            loss = np.repeat(loss, number_of_sectors, 1)

        return loss

if __name__ == '__main__':
    from matplotlib import pyplot as plt

    distance= np.arange(15000 , 400000, step=1)
    plt.figure(1)
    a = PropagationFreeSpace.get_loss(self="self", distance_2D= 15000, distance_3D = 0, frequency = 2550, number_of_sectors = 1)
    b = PropagationFreeSpace.get_loss(self="self", distance_2D= 20000, distance_3D = 0, frequency = 2550, number_of_sectors = 1)
    print(a)
    print(b)
    a_min1 = 15
    mass1 = 124.10
    a_min2 = 20
    mass2 = 126.60
    plt.plot(a_min1, mass1, 'r.', label='15 km / 124.10 dB')
    plt.plot(a_min2, mass2, 'o-', color = 'black', label='20 km / 126.10 dB')

    i = np.arange(15000, 400000, step=1)
    plt.legend(loc='upper right')
    plt.semilogx(i/1000, PropagationFreeSpace.get_loss(self="self", distance_2D= distance, distance_3D = 0, frequency = 2550, number_of_sectors = 1), 'b-', label='Free Space Path Loss')
    plt.grid(which='minor', alpha=0.2)
    plt.grid(which='major', alpha=0.5)
    plt.grid(True, color='b', linestyle='--', linewidth=0.2)
    plt.title('Free space path loss - Comparison (2D Distance)')
    plt.xlabel('Distância [km]')
    plt.ylabel('Atenuação [dB]')
    plt.legend()
    plt.show()
