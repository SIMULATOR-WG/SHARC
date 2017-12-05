# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 11:06:56 2017

@author: Calil
"""

from sharc.support.enumerations import StationType

import numpy as np
import matplotlib.pyplot as plt

class SpectralMaskImt(object):
    
    def __init__(self,sta_type: StationType, p_tx: float, freq_mhz: float, scenario = "OUTDOOR"):
        
        # Convert frequncy (for local use only)
        freq_ghz = freq_mhz/1e3
        
        # Spurious domain limits [dDm/MHz]
        self.spurious_limits = -13
        
        # Mask delta f breaking limits [MHz]
        self.delta_f_lim = np.array([0, 20, 400])
        
        if sta_type is StationType.IMT_UE:
            self.mask_dbm = np.array([p_tx, -5, -13, self.spurious_limits])
            
        elif sta_type is StationType.IMT_BS:
            if scenario is "OUTDOOR":
                
                # Table 4
                if p_tx > 32.5 and (freq_ghz > 37 and freq_ghz < 52.6):
                    self.mask_dbm = np.array([p_tx, -5, -13, self.spurious_limits])
                
                # Table 5
                elif p_tx < 32.5 and (freq_ghz > 37 and freq_ghz < 52.6):
                    self.mask_dbm = np.array([p_tx, -5, np.max((p_tx-45.5,-20)), 
                                              self.spurious_limits])
        
        # Spectral mask in dBc
        self.mask_dbc = self.mask_dbm - p_tx
        
        
    def power_calc(self,delta_f: float, band: float):
        
        # Limit delta f
        lim_df = delta_f + band
        
        # Included delta f values
        inc_df = np.where(np.logical_and(self.delta_f_lim > delta_f,
                                         self.delta_f_lim < lim_df))[0]
        
        # Power in mW
        if len(inc_df) == 0: 
            msk = self.mask_dbm[np.where(self.delta_f_lim > lim_df)]
            if len(msk) == 0: msk = np.array([self.mask_dbm[-1]])
            power = band*np.power(10,(msk[0]/10))
        else:
            power = (self.delta_f_lim[inc_df[0]] - delta_f)*\
                    np.power(10,(self.mask_dbm[inc_df[0]]/10)) + \
                    (lim_df - self.delta_f_lim[inc_df[-1]])*\
                    np.power(10,(self.mask_dbm[inc_df[-1] + 1]/10))
                
            for df in inc_df[0:-1]:
                power += (self.delta_f_lim[df + 1] - self.delta_f_lim[df])*\
                          np.power(10,(self.mask_dbm[df + 1]/10))
                    
                    
        return 10*np.log10(power)
        
            
        
if __name__ == '__main__':
    # Initialize variables
    sta_type = StationType.IMT_BS
    p_tx = 25.1
    freq = 43000
    
    # Create mask
    msk = SpectralMaskImt(sta_type,p_tx,freq)
    
    # Delta f
    delta_f = np.linspace(-100,600,num=1000)
    
    # Mask values
    mask_val = msk.mask_dbc[-1]*np.ones_like(delta_f)
    mask_val[np.where(delta_f < msk.delta_f_lim[-1])] = msk.mask_dbc[-2]
    mask_val[np.where(delta_f < msk.delta_f_lim[-2])] = msk.mask_dbc[-3]
    mask_val[np.where(delta_f < msk.delta_f_lim[-3])] = 0
    
    # Plot
    plt.plot(delta_f,mask_val)
    plt.xlim([delta_f[0],delta_f[-1]])
    plt.xlabel("$\Delta$f [MHz]")
    plt.ylabel("Spectral Mask [dBc]")
    plt.grid()
    plt.show()
            