# -*- coding: utf-8 -*-
"""
Created on Fri Mar  2 14:28:25 2018

@author: Calil
"""

import numpy as np
from itertools import product
from scipy.integrate import dblquad

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
        self.phi_vals = np.arange(-180,180,res)
        self.theta_vals = np.arange(0,180,res)
        self.antenna = None  
               
    def generate_correction_matrix(self, par: AntennaPar, c_chan: bool):
        """
        
        """
        # Create antenna object
        azi = 0 # Antenna azimuth: 0 degrees for simplicity
        ele = 0 # Antenna elevation: 0 degrees as well
        self.antenna = AntennaBeamformingImt(par,azi,ele)
        
        # Loop throug all the possible beams
        for phi, theta in product(self.phi_vals,self.theta_vals):
            pass
            
    def calculate_correction_factor(self, phi_e: float, theta_t: float, c_chan: bool):
        """
        
        """
        if c_chan:
            self.antenna.add_beam(phi_e,theta_t)
            beam = np.array([len(self.antenna.beams_list) - 1])
            int_f = lambda t,p: self.antenna._beam_gain(p,t,beam)*np.sin(t)
        else:
            int_f = lambda t,p: self.element.element_pattern(p,t)*np.sin(t)
        
        int_val = dblquad(int_f,-180,180,lambda t: 0, lambda t: 180)
        
        return int_val/(4*np.pi)
    
    def _save_files(self):
        pass
    
if __name__ == '__main__':
    pass