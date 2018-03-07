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
    """
    Implements spectral masks from document ITU 265-E. The masks are in the
    document's tables 1 to 8.
    
    Attributes:
        spurious_emissions (float): level of power emissions at spurious
            domain (dBm/MHz). Hardcoded as -13 dBm/MHz,  as specified in
            document ITU 265-E
        delta_f_lin (np.array): mask delta f breaking limits in MHz. Delta f 
            values for which the spectral mask changes value. Hard coded as
            [0, 20, 400]. In this context, delta f is the frequency distance to
            the transmission's edge frequencies
        freq_lim (no.array): frequency values for which the spectral mask
            changes emission value
        sta_type (StationType): type of station to which consider the spectral
            mask. Possible values are StationType.IMT_BS and StationType.IMT_UE
        freq_mhz (float): center frequency of station in MHz
        band_mhs (float): transmitting bandwidth of station in MHz
        scenario (str): INDOOR or OUTDOOR scenario
        p_tx (float): station's transmit power in dBm/MHz
        mask_dbm (np.array): spectral mask emission values in dBm
    """
    def __init__(self,sta_type: StationType, freq_mhz: float, band_mhz: float, scenario = "OUTDOOR"):
        """
        Class constructor.
        
        Parameters:
            sta_type (StationType): type of station to which consider the spectral
                mask. Possible values are StationType.IMT_BS and StationType.
                IMT_UE
            freq_mhz (float): center frequency of station in MHz
            band_mhs (float): transmitting bandwidth of station in MHz
            scenario (str): INDOOR or OUTDOOR scenario
        """
        # Spurious domain limits [dDm/MHz]
        self.spurious_emissions = -13
        # Mask delta f breaking limits [MHz]
        self.delta_f_lim = np.array([0, 20, 400])
        
        # Attributes
        self.sta_type = sta_type
        self.scenario = scenario
        self.band_mhz = band_mhz
        self.freq_mhz = freq_mhz
        
        self.freq_lim = np.concatenate(((freq_mhz - band_mhz/2)-self.delta_f_lim[::-1],
                                        (freq_mhz + band_mhz/2)+self.delta_f_lim))
       
        
    def set_mask(self, power = 0):
        """
        Sets the spectral mask (mask_dbm attribute) based on station type, 
        operating frequency and transmit power.
        
        Parameters:
            power (float): station transmit power. Default = 0
        """
        self.p_tx = power - 10*np.log10(self.band_mhz)
        
        # Set new transmit power value       
        if self.sta_type is StationType.IMT_UE:
            # Table 8
            mask_dbm = np.array([-5, -13, self.spurious_emissions])
            
        elif self.sta_type is StationType.IMT_BS and self.scenario is "INDOOR":             
            # Table 1
            mask_dbm = np.array([-5, -13, self.spurious_emissions])
            
        else:
            
            if (self.freq_mhz > 24250 and self.freq_mhz < 33400):
                if power >= 34.5:
                    # Table 2
                    mask_dbm = np.array([-5, -13, self.spurious_emissions])
                else:
                    # Table 3
                    mask_dbm = np.array([-5, np.max((power-47.5,-20)), 
                                          self.spurious_emissions])
            elif (self.freq_mhz > 37000 and self.freq_mhz < 52600):
                if power >= 32.5:
                    # Table 4
                    mask_dbm = np.array([-5, -13, self.spurious_emissions])
                else:
                    # Table 5
                    mask_dbm = np.array([-5, np.max((power-45.5,-20)), 
                                          self.spurious_emissions])
            elif (self.freq_mhz > 66000 and self.freq_mhz < 86000):
                if power >= 30.5:
                    # Table 6
                    mask_dbm = np.array([-5, -13, self.spurious_emissions])
                else:
                    # Table 7
                    mask_dbm = np.array([-5, np.max((power-43.5,-20)), 
                                          self.spurious_emissions])
            else:
                # Dummy spectral mask, for testing purposes only
                mask_dbm = np.array([-10, -20, -50])
                 
        self.mask_dbm = np.concatenate((mask_dbm[::-1],np.array([self.p_tx]),
                                        mask_dbm))
        
if __name__ == '__main__':
    # Initialize variables
    sta_type = StationType.IMT_BS
    p_tx = 25.1
    freq = 43000
    band = 200
    
    # Create mask
    msk = SpectralMaskImt(sta_type,freq,band)
    msk.set_mask(p_tx)
    
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
