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

