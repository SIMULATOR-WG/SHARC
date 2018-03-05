# -*- coding: utf-8 -*-
"""
Created on Fri Mar  2 14:28:25 2018

@author: Calil
"""

import numpy as np
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
        self.phi_min = -180
        self.phi_max = 180
        self.theta_min = 0
        self.theta_max = 180
        self.phi_vals = np.arange(-180,180,res)
        self.theta_vals = np.arange(0,180,res)
        self.antenna = None  
               
    def generate_correction_matrix(self, par: AntennaPar, c_chan: bool, file_name: str):
        """
        
        """
        # Create antenna object
        azi = 0 # Antenna azimuth: 0 degrees for simplicity
        ele = 0 # Antenna elevation: 0 degrees as well
        self.antenna = AntennaBeamformingImt(par,azi,ele)
        
        if c_chan:
            # Correction factor numpy array
            cf = np.zeros((len(self.phi_vals),len(self.theta_vals)))
            
            # Loop throug all the possible beams
            for phi_idx, phi in enumerate(self.phi_vals):
                for theta_idx, theta in enumerate(self.theta_vals):
                    cf[phi_idx,theta_idx] = self.calculate_correction_factor(phi,theta,c_chan)
                
        else:
            # Correction factor float
            cf = self.calculate_correction_factor(0,0,c_chan)
          
        # Convert to dB and save
        cf = 10*np.log10(cf)
        self._save_files(cf,par,file_name)
            
    def calculate_correction_factor(self, phi_e: float, theta_t: float, c_chan: bool):
        """
        
        """
        if c_chan:
            self.antenna.add_beam(phi_e,theta_t)
            beam = np.array([len(self.antenna.beams_list) - 1])
            int_f = lambda t,p: np.power(10,self.antenna._beam_gain(p,t,beam)/10)*np.sin(t)
        else:
            int_f = lambda t,p: np.power(10,self.element.element_pattern(p,t)/10)*np.sin(t)
        
        int_val = dblquad(int_f,self.phi_min,self.phi_max,
                          lambda t: self.theta_min, lambda t: self.theta_max)
        
        return int_val/(4*np.pi)
    
    def _save_files(self, correction, par:AntennaPar, file_name: str):
        """
        
        """
        np.savez(file_name,
                 correction_factor = correction,
                 parameters = np.array(par))
    
if __name__ == '__main__':
    pass