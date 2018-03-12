# -*- coding: utf-8 -*-
"""
Created on Fri Mar  2 14:28:25 2018

@author: Calil
"""

import numpy as np
from scipy.integrate import dblquad

from sharc.antenna.antenna_beamforming_imt import AntennaBeamformingImt
from sharc.support.named_tuples import AntennaPar

 
class BeamformingNormalizer(object):
    """
    Calculates the total integrated gain of given antenna array or element and
    returns its correction factor.
    
    Attributes:
        res_deg (float): correction factor matrix resolution [degreees]
        tol (float): integral absolute tolerance
        pni_min_deg (float): minimum phi escan value for correction factor
            matrix [degrees]. Hardcoded as -180 degrees
        phi_max_deg (float): maximum phi escan value for correction factor
            matrix [degrees]. Hardcoded as +180 degrees
        theta_min_deg (float): minimum theta tilt value for correction factor
            matrix [degrees]. Hardcoded as 0 degrees
        theta_max_deg (float): maximum theta tilt value for correction factor
            matrix [degrees]. Hardcoded as +180 degrees
        phi_min_rad (float): same as pni_min_deg, but in radians
        phi_max_rad (float): same as pni_max_deg, but in radians
        theta_min_rad (float): same as theta_min_deg, but in radians
        theta_max_rad (float): same as theta_max_deg, but in radians
        self.phi_vals_deg (np.array): all phi escan values for given resolution
        self.phi_vals_deg (np.array): all theta tilt values for given 
            resolution
        antenna (AntennaBeamformingImt): antenna to which calculate 
            normalization
    """
    def __init__(self, res_deg: float, tol: float):
        """
        Class constructor
        
        Parameters:
            res_deg (float): correction factor matrix resolution in degrees
            tol (float): absolute tolerance for integration
        """
        # Initialize attributes
        self.res_deg = res_deg
        self.tol = tol
        
        self.phi_min_deg = -180
        self.phi_max_deg = 180
        self.theta_min_deg = 0
        self.theta_max_deg = 180
        
        self.phi_min_rad = np.deg2rad(self.phi_min_deg)
        self.phi_max_rad = np.deg2rad(self.phi_max_deg)
        self.theta_min_rad = np.deg2rad(self.theta_min_deg)
        self.theta_max_rad = np.deg2rad(self.theta_max_deg)
        
        self.phi_vals_deg = np.arange(self.phi_min_deg,
                                      self.phi_max_deg,res_deg)
        self.theta_vals_deg = np.arange(self.theta_min_deg,
                                        self.theta_max_deg,res_deg)
        self.antenna = None  
               
    def generate_correction_matrix(self, par: AntennaPar, c_chan: bool, file_name: str):
        """
        Generates the correction factor matrix and saves it in a file
        
        Parameters:
            par (AntennaPar): set of parameters antenna parameters to be used
            c_chan (bool): if True, whole antenna array will be used. If false,
                just a single element is used
            file_name (str): name of file to which save the correction matrix
            
        Returns:
            cf (np.array): corraction factor matrix
            err (np.array): error matrix
        """
        # Create antenna object
        azi = 0 # Antenna azimuth: 0 degrees for simplicity
        ele = 0 # Antenna elevation: 0 degrees as well
        self.antenna = AntennaBeamformingImt(par,azi,ele)
        
        if c_chan:
            # Correction factor numpy array
            cf = np.zeros((len(self.phi_vals_deg),len(self.theta_vals_deg)))
            err = np.empty((len(self.phi_vals_deg),len(self.theta_vals_deg)), dtype=tuple)
            
            # Loop throug all the possible beams
            for phi_idx, phi in enumerate(self.phi_vals_deg):
                for theta_idx, theta in enumerate(self.theta_vals_deg):
                    cf[phi_idx,theta_idx], err[phi_idx,theta_idx] = self.calculate_correction_factor(phi,theta,c_chan)
                
        else:
            # Correction factor float
            cf, err = self.calculate_correction_factor(0,0,c_chan)
          
        # Save in file
        self._save_files(cf,par,file_name)
        return cf, err
            
    def calculate_correction_factor(self, phi_e: float, theta_t: float, c_chan: bool):
        """
        Calculates single correction factor
        
        Parameters:
            phi_e (float): phi escan value
            theta_t (float): theta tilt value
            c_chan (bool): if True, whole antenna array will be used. If false,
                just a single element is used
                
        Returns:
            cf (float): correction factor value
            err (tuple): upper and lower error bounds
        """
        if c_chan:
            self.antenna.add_beam(phi_e,theta_t)
            beam = int(len(self.antenna.beams_list) - 1)
            int_f = lambda t,p: \
            np.power(10,self.antenna._beam_gain(np.rad2deg(p),np.rad2deg(t),beam)/10)*np.sin(t)
        else:
            int_f = lambda t,p: \
            np.power(10,self.antenna.element.element_pattern(np.rad2deg(p),np.rad2deg(t))/10)*np.sin(t)
        
        int_val, err = dblquad(int_f,self.phi_min_rad,self.phi_max_rad,
                          lambda p: self.theta_min_rad, 
                          lambda p: self.theta_max_rad,
                          epsabs=self.tol,
                          epsrel=0.0)
        
        cf = -10*np.log10(int_val/(4*np.pi))
        
        hig_bound = -10*np.log10((int_val - err)/(4*np.pi))
        low_bound = -10*np.log10((int_val + err)/(4*np.pi))
        
        return cf, (low_bound,hig_bound)
    
    def _save_files(self, correction, par:AntennaPar, file_name: str):
        np.savez(file_name,
                 correction_factor = correction,
                 parameters = par)
    
if __name__ == '__main__':
    """
    Plots correction factor for horizontal and vertical planes.
    """
    import matplotlib.pyplot as plt
    
    # Create normalizer object
    resolution = 2.5
    tolerance = 1e-1
    norm = BeamformingNormalizer(resolution,tolerance)
    
    # Antenna parameters
    element_pattern = "M2101"
    element_max_g = 5
    element_phi_deg_3db = 65
    element_theta_deg_3db = 65
    element_am = 30
    element_sla_v = 30
    n_rows = 8
    n_columns = 8
    horiz_spacing = 0.5
    vert_spacing = 0.5
    down_tilt = 0
    par = AntennaPar(element_pattern,
                     element_max_g,
                     element_phi_deg_3db,
                     element_theta_deg_3db,
                     element_am,
                     element_sla_v,
                     n_rows,
                     n_columns,
                     horiz_spacing,
                     vert_spacing,
                     down_tilt)
    
    # Set range of values & calculate correction factor
    norm.theta_vals_deg = np.array([90])
    c_chan = True
    file_name = 'main_test.npz'
    cf, err = norm.generate_correction_matrix(par,c_chan,file_name)
    err_low, err_high = zip(*np.ravel(err))
    
    plt.plot(norm.phi_vals_deg,cf,
             norm.phi_vals_deg,err_low,'r--',
             norm.phi_vals_deg,err_high,'r--')
    plt.xlim(-180,180)
    plt.ylabel(r"Correction factor [dB]")
    plt.xlabel(r"Azimuth angle $\phi$ [deg]")
    plt.title(r"Elevation angle $\theta$ = 90 deg")
    plt.show()
    
    # Set range of values & calculate correction factor
    norm = BeamformingNormalizer(resolution,tolerance)
    norm.phi_vals_deg = np.array([0])
    c_chan = True
    file_name = 'main_test.npz'
    cf, err = norm.generate_correction_matrix(par,c_chan,file_name)
    err_low, err_high = zip(*np.ravel(err))
    
    plt.plot(norm.theta_vals_deg,np.transpose(cf),
             norm.theta_vals_deg,np.transpose(err_low),'r--',
             norm.theta_vals_deg,np.transpose(err_high),'r--')
    plt.xlim(0,180)
    plt.ylabel(r"Correction factor [dB]")
    plt.xlabel(r"Elevation angle $\theta$ [deg]")
    plt.title(r"Azimuth angle $\phi$ = 0 deg")
    plt.show()
    
    