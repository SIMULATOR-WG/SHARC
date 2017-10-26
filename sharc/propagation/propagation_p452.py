# -*- coding: utf-8 -*-
"""
Created on Mon Apr 17 15:35:00 2017

@author: andre barreto
"""

from sharc.propagation.propagation import Propagation

import math
import numpy as np
from sharc.propagation.propagation_free_space import PropagationFreeSpace
from sharc.propagation.propagation_clutter_loss import PropagationClutterLoss
from sharc.propagation.propagation_building_entry_loss import PropagationBuildingEntryLoss
from sharc.propagation.atmosphere import ReferenceAtmosphere
from sharc.parameters.parameters_fss_ss import ParametersFssSs
from sharc.support.enumerations import StationType
from sharc.propagation.scintillation import Scintillation

class Propagation452(Propagation):
    """
    Implements the earth-to-space channel model from ITU-R P.619

    Public methods:
        get_loss: Calculates path loss for earth-space link
    """

    def __init__(self):
        super().__init__()


    def get_loss(self, *args, **kwargs) -> np.array:
        """
        Calculates path loss for earth-space link

        Parameters
        ----------
            distance_3D (np.array) : distances between stations [m]
            frequency_MHz (np.array) : center frequencies [MHz]
            indoor_stations (np.array) : array indicating stations that are indoor
            percentage_loss (float/string) : required  time  percentage  for  which  the  calculated  basic
                                             transmission loss is not exceeded
            tx_lat, tx_long (np.array) : latitude and longitude of transmitter
            rx_lat, rx_long (np.array) : latitude and longitude of receiver


        Returns
        -------
            array with path loss values with dimensions of distance_3D

        """

        distance_3D = kwargs["distance_3D"]
        frequency_MHz = kwargs["frequency_MHz"]
        freq_GHz = frequency_MHz / 1000
        indoor_stations = kwargs["indoor_stations"]

        loss = 0
        return loss

if __name__ == '__main__':
    a = "dummy"
