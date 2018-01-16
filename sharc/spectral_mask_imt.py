# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 11:06:56 2017

@author: Calil
"""

from sharc.support.enumerations import StationType

import numpy as np
import matplotlib.pyplot as plt

class SpectralMaskImt(object):
    
    def __init__(self,sta_type: StationType, freq_mhz: float, band_mhz: float, scenario = "OUTDOOR"):
        
        # Spurious domain limits [dDm/MHz]
        self.spurious_limits = -13
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
        
        self.p_tx = power - 10*np.log10(self.band_mhz)
        
        # Set new transmit power value       
        if self.sta_type is StationType.IMT_UE:
            # Table 8
            mask_dbm = np.array([-5, -13, self.spurious_limits])
            
        elif self.sta_type is StationType.IMT_BS and self.scenario is "INDOOR":             
            # Table 1
            mask_dbm = np.array([-5, -13, self.spurious_limits])
            
        else:
            
            if (self.freq_mhz > 24250 and self.freq_mhz < 33400):
                if self.p_tx >= 34.5:
                    # Table 2
                    mask_dbm = np.array([-5, -13, self.spurious_limits])
                else:
                    # Table 3
                    mask_dbm = np.array([-5, np.max((power-47.5,-20)), 
                                          self.spurious_limits])
            elif (self.freq_mhz > 37000 and self.freq_mhz < 52600):
                if self.p_tx >= 32.5:
                    # Table 4
                    mask_dbm = np.array([-5, -13, self.spurious_limits])
                else:
                    # Table 5
                    mask_dbm = np.array([-5, np.max((power-45.5,-20)), 
                                          self.spurious_limits])
            elif (self.freq_mhz > 66000 and self.freq_mhz < 86000):
                if self.p_tx >= 30.5:
                    # Table 6
                    mask_dbm = np.array([-5, -13, self.spurious_limits])
                else:
                    # Table 7
                    mask_dbm = np.array([-5, np.max((power-43.5,-20)), 
                                          self.spurious_limits])
            else:
                # Dummy spectral mask, for testing purposes only
                mask_dbm = np.array([-10, -20, -50])
                 
        self.mask_dbm = np.concatenate((mask_dbm[::-1],np.array([self.p_tx]),
                                        mask_dbm))
        
    def power_calc(self,center_f: float, band: float):
        
        # Limit delta f
        df_min = center_f - band/2
        df_max = center_f + band/2
        
        # Power in mW
        power_oob = 0 # Out-of-band power
        
        # Included delta f values
        inc_df = np.where(np.logical_and(self.freq_lim > df_min,
                                         self.freq_lim < df_max))[0]
        if len(inc_df) == 0:
            
            msk = self.mask_dbm[np.where(self.freq_lim >= df_max)]
            if len(msk) == 0: msk = np.array([self.mask_dbm[-1]])
            pwr_lvl = msk[0]
            
            if pwr_lvl != self.p_tx: power_oob += band*np.power(10,(pwr_lvl/10))
            
        else:
            
            pwr_lvl_1 = self.mask_dbm[inc_df[0]]
            pwr_lvl_2 = self.mask_dbm[inc_df[-1] + 1]
                    
            if pwr_lvl_1 != self.p_tx: 
                power_oob += (self.freq_lim[inc_df[0]] - df_min)*\
                    np.power(10,(pwr_lvl_1/10))
                    
            if pwr_lvl_2 != self.p_tx: 
                power_oob += (df_max - self.freq_lim[inc_df[-1]])*\
                    np.power(10,(pwr_lvl_2/10))
                
            for df in inc_df[0:-1]:
                
                pwr_lvl = self.mask_dbm[df + 1]
                
                if pwr_lvl != self.p_tx: 
                    power_oob += (self.freq_lim[df + 1] - self.freq_lim[df])*\
                          np.power(10,(pwr_lvl/10))
                    
                    
        return 10*np.log10(power_oob)
        
            
        
if __name__ == '__main__':
    # Initialize variables
    sta_type = StationType.IMT_BS
    p_tx = 25.1
    freq = 43000
    band = 200
    
    # Create mask
    msk = SpectralMaskImt(sta_type,freq,band)
    msk.set_power(p_tx)
    
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
            