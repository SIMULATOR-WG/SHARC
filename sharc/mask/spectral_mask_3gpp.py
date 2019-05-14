# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 11:06:56 2017

@author: Calil
"""

from sharc.support.enumerations import StationType
from sharc.mask.spectral_mask import SpectralMask

import numpy as np

class SpectralMask3Gpp(SpectralMask):

    def __init__(self,
                 sta_type: StationType, 
                 freq_mhz: float, 
                 band_mhz: float, 
                 spurious_emissions: float,
                 scenario = "OUTDOOR"):
        """
        Implements spectral mask from 3GPP 36.104 Table 6.6.3.1-6, which 
        defines requirements for Wide Area BS operating with 5, 10, 15 or 
        20 MHz channel bandwidth.
        In order to characterize Cat-A or Cat-B devices it is necessary to 
        adjust the input parameter that defines the spurious emission level:
            Cat-A: -13 dBm/MHz
            Cat-B: -30 dBm/MHz
        """
        # Spurious domain limits [dBm/MHz]
        self.spurious_emissions = spurious_emissions
        # Mask delta f breaking limits [MHz]
        self.delta_f_lim = np.arange(0, 5.1, .1)
        self.delta_f_lim = np.append(self.delta_f_lim, 10)

        # Attributes
        self.sta_type = sta_type
        self.scenario = scenario
        self.band_mhz = band_mhz
        self.freq_mhz = freq_mhz

        #delta_f_lim_flipped = np.flip(self.delta_f_lim,0)
        delta_f_lim_flipped = self.delta_f_lim[::-1]
        
        self.freq_lim = np.concatenate(((freq_mhz - band_mhz/2) - delta_f_lim_flipped,
                                        (freq_mhz + band_mhz/2) + self.delta_f_lim))

    def set_mask(self, power = 0):
        self.p_tx = power - 10 * np.log10(self.band_mhz)
        # mask in dBm/MHz
        mask_dbm = 3 - 7 / 5 * (np.arange(.05, 5, .1) - .05)
        mask_dbm = np.append(mask_dbm, np.array([-4, self.spurious_emissions]))
        
        #mask_dBm_flipped = np.flip(mask_dbm, 0)
        mask_dBm_flipped = mask_dbm[::-1]

        self.mask_dbm = np.concatenate((mask_dBm_flipped, 
                                        np.array([self.p_tx]),
                                        mask_dbm))

