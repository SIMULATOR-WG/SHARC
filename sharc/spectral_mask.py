# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 11:06:56 2017

@author: Calil
"""

from abc import ABC, abstractmethod
import numpy as np

class SpectralMask(ABC):

    @abstractmethod
    def set_mask(self, p_tx = 0):
        pass

    def power_calc(self,center_f: float, band: float):
        """
        Calculates out-of-band power in the given band. It does that by 
        dividing the band into the rectangular sections defined by the spectral
        mask and adding up the area of all the rectangles.
        
        Parameters:
            center_f (float): center frequency of band in which out-of-band
                power is to be calculated
            band (float): bandwidth of band in which out-of-band power is to
                be calculated
        """
        # Limit delta f: edges of band
        df_min = center_f - band/2
        df_max = center_f + band/2

        # Power in mW
        power_oob = 0 # Out-of-band power

        # Included delta f values: values of spectral mask delta f break limist
        # which are contained within the band
        inc_df = np.where(np.logical_and(self.freq_lim > df_min,
                                         self.freq_lim < df_max))[0]
        
        # If no break limits are within band: the band does not need to be
        # divided into rectangles
        if len(inc_df) == 0:

            # Define what is the power emission level at that frequency
            msk = self.mask_dbm[np.where(self.freq_lim >= df_max)]
            # If df_max is below smallest break limit
            if len(msk) == 0: msk = np.array([self.mask_dbm[-1]])
            # Turn it into scalar
            pwr_lvl = msk[0]

            if pwr_lvl != self.p_tx: power_oob += band*np.power(10,(pwr_lvl/10))
        # If one or more break limitas are within band
        else:
            
            # Lower and upper power emission levels
            pwr_lvl_1 = self.mask_dbm[inc_df[0]]
            pwr_lvl_2 = self.mask_dbm[inc_df[-1] + 1]

            # Upper and lower rectangles
            if pwr_lvl_1 != self.p_tx:
                power_oob += (self.freq_lim[inc_df[0]] - df_min)*\
                    np.power(10,(pwr_lvl_1/10))
            if pwr_lvl_2 != self.p_tx:
                power_oob += (df_max - self.freq_lim[inc_df[-1]])*\
                    np.power(10,(pwr_lvl_2/10))

            # Middle rectangles
            for df in inc_df[0:-1]:
                pwr_lvl = self.mask_dbm[df + 1]

                if pwr_lvl != self.p_tx:
                    power_oob += (self.freq_lim[df + 1] - self.freq_lim[df])*\
                          np.power(10,(pwr_lvl/10))

        return 10*np.log10(power_oob)

