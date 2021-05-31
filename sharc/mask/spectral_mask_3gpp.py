# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 11:06:56 2017

@author: Calil

@modified: Luciano Camilo Tue Jan 26 13:49:25 2021
"""

from sharc.support.enumerations import StationType
from sharc.mask.spectral_mask import SpectralMask

import numpy as np
import matplotlib.pyplot as plt
import sys


class SpectralMask3Gpp(SpectralMask):

    def __init__(self,
                 sta_type: StationType,
                 freq_mhz: float,
                 band_mhz: float,
                 spurious_emissions: float,
                 scenario="HIBS"):
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
        # Spurious domain limits [dDm/MHz]
        self.spurious_limits = -13
        # Mask delta f breaking limits [MHz]
        self.delta_f_lim = np.arange(0, 5.1, .1)
        self.delta_f_lim = np.append(self.delta_f_lim, 10)

        if sta_type is not StationType.IMT_BS and sta_type is not StationType.IMT_UE:
            message = "ERROR\nInvalid station type: " + str(sta_type)
            sys.stderr.write(message)
            sys.exit(1)

        if band_mhz not in [5, 10, 15, 20]:
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
        # delta_f_lim_flipped = np.flip(self.delta_f_lim,0)
        delta_f_lim_flipped = delta_f_lim[::-1]

        self.freq_lim = np.concatenate(((self.freq_mhz - self.band_mhz/2) - delta_f_lim_flipped,
                                        (self.freq_mhz + self.band_mhz/2) + delta_f_lim))

    def get_frequency_limits(self,
                             sta_type: StationType,
                             bandwidth: float) -> np.array:
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

    def set_mask(self, power=0):
        emission_limits = self.get_emission_limits(self.sta_type,
                                                   self.band_mhz,
                                                   self.spurious_emissions)
        # print(emission_limits)
        self.p_tx = power - 10 * np.log10(self.band_mhz)
        emission_limits_flipped = emission_limits[::-1]
        self.mask_dbm = np.concatenate((emission_limits_flipped,
                                        np.array([self.p_tx]),
                                        emission_limits))

    def get_emission_limits(self,
                            sta_type: StationType,
                            bandwidth: float,
                            spurious_emissions: float) -> np.array:
        if sta_type is StationType.IMT_BS:
            # emission limits in dBm/MHz
            emission_limits = 3 - 7 / 5 * (np.arange(.05, 5, .1) - .05)
            emission_limits = np.append(emission_limits, np.array([-4, spurious_emissions]))
            # emission_limits = np.append(emission_limits, np.array([-4, -60]))
            # print(emission_limits)
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


if __name__ == '__main__':
    # Initialize variables
    sta_type = StationType.IMT_BS
    p_tx = 43
    freq = 2680 # 5 MHz guard band RAS
    band = 20

    # Create mask
    msk = SpectralMask3Gpp(sta_type, freq, band, -13)
    msk.set_mask(p_tx)

    # Frequencies
    freqs = np.linspace(-40, 40, num=1000)+freq

    # Mask values
    mask_val = np.ones_like(freqs)*msk.mask_dbm[0]
    for k in range(len(msk.freq_lim)-1, -1, -1):
        mask_val[np.where(freqs < msk.freq_lim[k])] = msk.mask_dbm[k]

    # Plot

    rasx = np.linspace(2700, 2701.5, 60)
    rasy = np.linspace(0, 0, 60)
    rasx_ = np.linspace(2700, 2700, 60)
    rasy_ = np.linspace(-13, 0, 60)
    rasx_1 = np.linspace(2701.5, 2701.5, 60)
    rasy_1 = np.linspace(0, -13, 60)

    plt.plot(rasx_1, rasy_1, 'r--', linewidth=2, color = 'orange', label= ' Meteorological Radar Band ')
    plt.plot(rasx_, rasy_, 'r--', linewidth=2, color='orange')
    plt.plot(rasx, rasy, 'r--', linewidth=2, color='orange')
    plt.plot(freqs, mask_val, 'r-', linewidth=1.5, color='black', label=' HIBS Spectral Mask ')
    plt.legend(loc='upper right')
    # plt.grid(which='minor', alpha=0.7)
    # plt.grid(which='major', alpha=0.7)
    # plt.grid(True, color='k', linestyle='--', linewidth=0.5)
    plt.xlim([freqs[0], freqs[-1]])
    #plt.xlim([freqs[0], 2730])
    # plt.title('3GPP 36.104 Spectral Mask')
    plt.xlabel("$\Delta$f [MHz]")
    plt.ylabel("Spectral Mask [dBc]")
    # plt.grid()
    plt.show()
