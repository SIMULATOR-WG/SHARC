# -*- coding: utf-8 -*-
"""
Created on Fri Mar  2 14:28:25 2018

@author: Calil
"""

import numpy as np

from sharc.antenna.antenna_beamforming_imt import AntennaBeamformingImt
from sharc.support.named_tuples import AntennaPar

 
class BeamformingNormalization(object):
    """
    
    """
    def _init_(self, res: float):
        """
        
        """
        # Initialize attributes
        self.res = res
        self.phi_scan_vals = np.arange(-180,180,res)
        self.theta_tilt_vals = np.arange(0,180,res)
        self.antenna = None  
               
    def calculcate_correction_factor(self):
        """
        
        """
        pass
    
    def calculate_gain_matrix(self, par: AntennaPar, beams_list: np.array, c_chan: bool):
        """
        
        """
        # Create antenna object
        azi = 0 # Antenna azimuth: 0 degrees for simplicity
        ele = 0 # Antenna elevation: 0 degrees as well
        self.antenna = AntennaBeamformingImt(par,azi,ele)
        
        gains = np.zeros((len(self.phi_scan_vals),len(self.theta_tilt_vals)))
        
        for k,theta_tilt in enumerate(self.theta_tilt_vals):
            theta_v = theta_tilt*np.ones_like(self.phi_scan_vals)
            gains[k,:] = self.antenna.calcualte_gain(phi_vec = self.phi_scan_vals,
                                                     theta_vec=theta_v,
                                                     beams_l = beams_list,
                                                     co_channel = c_chan)   
        return gains
    
    def save_files(self):
        pass
    
if __name__ == '__main__':
    pass