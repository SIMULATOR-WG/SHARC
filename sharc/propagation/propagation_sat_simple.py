# -*- coding: utf-8 -*-
"""
Created on Thu Feb 16 12:04:27 2017

@author: edgar
"""

from sharc.propagation.propagation import Propagation
from sharc.propagation.propagation_free_space import PropagationFreeSpace
from sharc.propagation.propagation_clutter_loss import PropagationClutterLoss
from sharc.support.enumerations import StationType

import numpy as np

class PropagationSatSimple(Propagation):
    """
    Implements the simplified satellite propagation model
    """

    def __init__(self):
        super().__init__()
        self.clutter = PropagationClutterLoss()
        self.free_space = PropagationFreeSpace()
        self.atmospheric_loss = 1
        self.polarization_loss = 3
        self.building_loss = 20

    def get_loss(self, *args, **kwargs) -> np.array:
        d = kwargs["distance_3D"]
        f = kwargs["frequency"]
        p = kwargs["loc_percentage"]
        indoor_stations = kwargs["indoor_stations"]
        elevation = kwargs["elevation"]

        free_space_loss = self.free_space.get_loss(distance_3D=d,
                                                   frequency=f)
        clutter_loss = np.maximum(0, self.clutter.get_loss(frequency=f,
                                                           elevation=elevation["free_space"],
                                                           loc_percentage=p,
                                                           station_type=StationType.FSS_SS))
        building_loss = self.building_loss*indoor_stations

        loss = (free_space_loss + clutter_loss + building_loss +
                self.polarization_loss + self.atmospheric_loss)

        return loss
