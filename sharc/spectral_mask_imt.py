# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 11:06:56 2017

@author: Calil
"""

from sharc.support.enumerations import StationType
from sharc.spectral_mask import SpectralMask

import numpy as np
import matplotlib.pyplot as plt

class SpectralMaskImt(SpectralMask):

    def __init__(self,sta_type: StationType, p_tx: float, freq_mhz: float, band_mhz: float, scenario = "OUTDOOR"):

        # Spurious domain limits [dDm/MHz]
        self.spurious_limits = -13
        # Mask delta f breaking limits [MHz]
        self.delta_f_lim = np.array([0, 20, 400])

        # Attributes
        self.tx_power = p_tx
        self.sta_type = sta_type
        self.scenario = scenario
        self.band_mhz = band_mhz
        self.freq_mhz = freq_mhz

        self.freq_lim = np.concatenate(((freq_mhz - band_mhz/2)-np.flip(self.delta_f_lim,0),
                                        (freq_mhz + band_mhz/2)+self.delta_f_lim))

        self.set_power(p_tx)


    def set_power(self, p_tx: float):

        self.p_tx = p_tx - 10*np.log10(self.band_mhz)

        # Set new transmit power value
        if self.sta_type is StationType.IMT_UE:
            mask_dbm = np.array([-5, -13, self.spurious_limits])

        elif self.sta_type is StationType.IMT_BS and self.scenario is "OUTDOOR":
            # Table 4
            if self.p_tx > 32.5 and (self.freq_mhz > 37000 and self.freq_mhz < 52600):
                mask_dbm = np.array([-5, -13, self.spurious_limits])

            # Table 5
            elif self.p_tx < 32.5 and (self.freq_mhz > 37000 and self.freq_mhz < 52600):
                mask_dbm = np.array([-5, np.max((p_tx-45.5,-20)),
                                          self.spurious_limits])

            else: # Dummy spectral mask, for testing purposes only
                mask_dbm = np.array([-10, -20, -50])

        self.mask_dbm = np.concatenate((np.flip(mask_dbm,0),np.array([self.p_tx]),
                                        mask_dbm))



if __name__ == '__main__':
    # Initialize variables
    sta_type = StationType.IMT_BS
    p_tx = 25.1
    freq = 43000

    # Create mask
    msk = SpectralMaskImt(sta_type,p_tx,freq)

    # Frequencies
    freqs = np.linspace(-600,600,num=1000)+freq

    # Mask values
    mask_val = np.ones_like(freqs)*msk.mask_dbm[0]
    for k in range(len(msk.freq_lim)-1,-1,-1):
        mask_val[np.where(freqs < msk.freq_lim[k])] = msk.mask_dbm[k]

    # Plot
    plt.plot(freqs,mask_val)
    plt.xlim([freqs[0],freqs[-1]])
    plt.xlabel("$\Delta$f [MHz]")
    plt.ylabel("Spectral Mask [dBc]")
    plt.grid()
    plt.show()
