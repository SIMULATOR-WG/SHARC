# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 11:06:56 2017

@author: Calil
"""

from sharc.support.enumerations import StationType
from sharc.mask.spectral_mask import SpectralMask

import numpy as np
import sys

class SpectralMask3Gpp(SpectralMask):

    def __init__(self,
                 sta_type: StationType, 
                 freq_mhz: float, 
                 band_mhz: float, 
                 spurious_emissions: float,
                 scenario = "OUTDOOR"):
        """
        Implements spectral emission mask from 3GPP 36.104 Table 6.6.3.1-6 for 
        Wide Area BS operating with 5, 10, 15 or 20 MHz channel bandwidth.
        
        Also implements spectral emission mask from 3GPP 36.101 Table 6.6.2.1.1-1 
        for UE operating with 5, 10, 15 or 20 MHz channel bandwidth.
        
        In order to characterize Cat-A or Cat-B base stations, it is necessary 
        to adjust the input parameter that defines the spurious emission level:
            Cat-A: -13 dBm/MHz
            Cat-B: -30 dBm/MHz
        """
        if sta_type is not StationType.IMT_BS and sta_type is not StationType.IMT_UE:
            message = "ERROR\nInvalid station type: " + str(sta_type)
            sys.stderr.write(message)
            sys.exit(1)
            
        if band_mhz not in [ 5, 10, 15, 20 ]:
            message = "ERROR\nInvalid bandwidth for 3GPP mask: " + band_mhz
            sys.stderr.write(message)
            sys.exit(1)             

        # Attributes
        self.spurious_emissions = spurious_emissions
        self.sta_type = sta_type
        self.scenario = scenario
        self.band_mhz = band_mhz
        self.freq_mhz = freq_mhz
        
        delta_f_lim = self.get_frequency_limits(self.sta_type, self.band_mhz)
        #delta_f_lim_flipped = np.flip(self.delta_f_lim,0)
        delta_f_lim_flipped = delta_f_lim[::-1]
        
        self.freq_lim = np.concatenate(((self.freq_mhz - self.band_mhz/2) - delta_f_lim_flipped,
                                        (self.freq_mhz + self.band_mhz/2) + delta_f_lim))

            
    def get_frequency_limits(self, 
                             sta_type : StationType,
                             bandwidth : float) -> np.array:
        """
        Calculates the frequency limits of the spectrum emission masks. This
        implementation is valid only for bandwidths equal to 5, 10, 15 or 20 MHz.
        """
        
        if sta_type is StationType.IMT_BS:
            # Mask delta f breaking limits [MHz]
            delta_f_lim = np.arange(0, 5.1, .1)
            delta_f_lim = np.append(delta_f_lim, 10)
        else:        
            delta_f_lim = np.array([0, 1, 5])
            if bandwidth == 5:
                delta_f_lim = np.append(delta_f_lim, np.array([6, 10]))
            if bandwidth == 10:
                delta_f_lim = np.append(delta_f_lim, np.array([10, 15]))
            if bandwidth == 15:
                delta_f_lim = np.append(delta_f_lim, np.array([15, 20]))
            else:
                delta_f_lim = np.append(delta_f_lim, np.array([20, 25]))
        return delta_f_lim
        
            
    def set_mask(self, power = 0):
        emission_limits = self.get_emission_limits(self.sta_type,
                                                   self.band_mhz,
                                                   self.spurious_emissions)
        self.p_tx = power - 10 * np.log10(self.band_mhz)
        #emission_limits = np.flip(emission_limits, 0)
        emission_limits_flipped = emission_limits[::-1]
        self.mask_dbm = np.concatenate((emission_limits_flipped, 
                                        np.array([self.p_tx]),
                                        emission_limits))


    def get_emission_limits(self, 
                            sta_type : StationType,
                            bandwidth : float,
                            spurious_emissions : float) -> np.array:
        if sta_type is StationType.IMT_BS:
            # emission limits in dBm/MHz
            emission_limits = 3 - 7 / 5 * (np.arange(.05, 5, .1) - .05)
            emission_limits = np.append(emission_limits, np.array([-4, spurious_emissions]))
        else:
            if bandwidth == 5:
                limit_r1 = np.array([-15])
            elif bandwidth == 10:
                limit_r1 = np.array([-18])
            elif bandwidth == 15:
                limit_r1 = np.array([-20])
            else:
                limit_r1 = np.array([-21])
            emission_limits = np.append(limit_r1 + 10*np.log10(1/0.03),
                                        np.array([-10, -13, -25, spurious_emissions]))
        return emission_limits