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
        p = kwargs["loc_percentage"]
        indoor_stations = kwargs["indoor_stations"]
        
        free_space_loss = self.free_space.get_loss(distance_2D = d,
                                                   frequency = f)
        clutter_loss_1 = np.zeros(d.shape)
        clutter_loss_2 = np.zeros(d.shape)
        
        # TODO: this implementation assumes that there is only one station of 
        # the other system and n IMT stations, so dimensions of the matrices
        # will always be 1 x n. This has to be generalized
        
        # when 250 < d < 1000, correction is applied at only one end of the path
        id_1 = np.where(d > 250)[1]
        if len(id_1):
            clutter_loss_1[0,id_1] = self.clutter.get_loss(frequency = f[0,id_1],
                                                         distance_2D = d[0,id_1],
                                                         loc_percentage = p,
                                                         station_type = StationType.FSS_ES)
            
        # when d > 1000, correction is applied at both ends of the path
        id_2 = np.where(d > 1000)[1]
        if len(id_2):
            clutter_loss_2[0,id_2] = self.clutter.get_loss(frequency = f[0,id_2],
                                                         distance_2D = d[0,id_2],
                                                         loc_percentage = p,
                                                         station_type = StationType.FSS_ES)
        
        building_loss = self.building_loss*indoor_stations
        
        loss = free_space_loss + building_loss + clutter_loss_1 + clutter_loss_2
        
        return loss
        

#if __name__ == '__main__':
#
#    import matplotlib.pyplot as plt
#        
#    ###########################################################################
#    # Print path loss for TerrestrialSimple and Free Space
#    
#    shadowing_std = 0
#    d = np.transpose(np.linspace(10, 10000, num=10000)[:,np.newaxis])
#    freq = 27000*np.ones(d.shape)
#    indoor_stations = np.zeros(d.shape, dtype = bool)
#    
#    free_space = PropagationFreeSpace()
#    ter_simple = PropagationTerSimple()
#    
#    loss_ter = ter_simple.get_loss(distance_2D = d, 
#                                  frequency = freq,
#                                  loc_percentage = 0.5,
#                                  indoor_stations = indoor_stations)
#
#    loss_fs = free_space.get_loss(distance_2D = d, 
#                                  frequency = freq)
#    
#    fig = plt.figure(figsize=(8,6), facecolor='w', edgecolor='k')
#
#    plt.semilogx(np.squeeze(d), np.squeeze(loss_fs), label = "free space")
#    plt.semilogx(np.squeeze(d), np.squeeze(loss_ter), label = "free space + clutter loss")
#    
#    plt.title("Free space with additional median clutter loss ($f=27GHz$)")
#    plt.xlabel("distance [m]")
#    plt.ylabel("path loss [dB]")
#    plt.xlim((0, d[-1,0]))
#    plt.ylim((80, 240))                
#    plt.legend(loc="upper left")
#    plt.tight_layout()    
#    plt.grid()
#
#    plt.show()
            