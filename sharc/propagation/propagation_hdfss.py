# -*- coding: utf-8 -*-
"""
Created on Mon Jul 30 17:28:47 2018

@author: Calil
"""

import numpy as np

from sharc.propagation.propagation import Propagation
from sharc.propagation.propagation_p1411 import PropagationP1411
from sharc.propagation.propagation_free_space import PropagationFreeSpace

class PropagationHDFSS(Propagation):
    """
    This is a wrapper class which can be used for indoor simulations. It
    calculates the basic path loss between IMT stations and HDFSS Earth Stations.
    It uses:
        FSPL for distance < 55m
        P.1411 LOS for distances 55m < distance < 260m
        P.1411 NLOS for distances distance > 260m
    """
    def __init__(self, random_number_gen: np.random.RandomState):
        super().__init__(random_number_gen)
        
        self.propagation_p1411 = PropagationP1411(random_number_gen)
        self.propagation_fspl = PropagationFreeSpace(random_number_gen)
        
        self.fspl_dist = 55
        self.p1411_los_dist = 260
        
    def get_loss(self, *args, **kwargs) -> np.array:
        """
        Calculates path loss for given distances and frequencies

        Parameters
        ----------
            distance_3D (np.array) : 3D distances between stations
            distance_2D (np.array) : 2D distances between stations
            frequency (np.array) : center frequencie [MHz]
            indoor (np.array) : indicates whether UE is indoor
            shadowing (bool) : if shadowing should be added or not
            number_of_sectors (int): number of sectors in a node (default = 1)

        Returns
        -------
            array with path loss values with dimensions of distance_2D

        """
        if "distance_3D" in kwargs:
            d = kwargs["distance_3D"]
        else:
            d = kwargs["distance_2D"]

        f = kwargs["frequency"]/1e3
        shad = kwargs.pop("shadow",True)
        number_of_sectors = kwargs.pop("number_of_sectors",1)
        
        fspl_idx = np.where(d < self.fspl_dist)
        p1411_los_idx = np.where(np.logical_and(d > self.fspl_dist,d < self.p1411_los_dist))
        p1411_nlos_idx = np.where(d > self.p1411_los_dist)
        
        loss = np.zeros_like(d)
        
        loss[fspl_idx] = self.propagation_fspl.get_loss(distance_3D=d[fspl_idx],
                                                        frequency=f[fspl_idx])
        loss[p1411_los_idx] = self.propagation_p1411.get_loss(distance_3D=d[p1411_los_idx],
                                                              frequency=f[p1411_los_idx],
                                                              los=True,
                                                              shadow=shad)
        loss[p1411_nlos_idx] = self.propagation_p1411.get_loss(distance_3D=d[p1411_nlos_idx],
                                                               frequency=f[p1411_nlos_idx],
                                                               los=False,
                                                               shadow=shad)
        
        if number_of_sectors > 1:
            loss = np.repeat(loss, number_of_sectors, 1)
            
        return loss
    
if __name__ == '__main__':
    
    import matplotlib.pyplot as plt
    
    rnd = np.random.RandomState(101)
    prop = PropagationHDFSS(rnd)
    
    d = np.linspace(5,1000,num=1000)
    f = 40e3*np.ones_like(d)
        
    loss = prop.get_loss(distance_3D=d,
                         frequency=f,
                         shadow=False)
    
    plt.plot(d,loss)
    plt.xlabel("Distance [m]")
    plt.ylabel("Path Loss [dB]")
    plt.grid()
    plt.show()