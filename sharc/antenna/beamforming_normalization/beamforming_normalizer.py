# -*- coding: utf-8 -*-
"""
Created on Fri Mar  2 14:28:25 2018

@author: Calil
"""

import numpy as np
from scipy.integrate import dblquad
from sys import stdout

from sharc.antenna.antenna_beamforming_imt import AntennaBeamformingImt
from sharc.support.named_tuples import AntennaPar


class BeamformingNormalizer(object):
    """
    Calculates the total integrated gain of given antenna array or element and
    returns its correction factor.

    Attributes:
        resolution_deg (float): correction factor matrix resolution [degreees]
        tolerance (float): integral absolute tolerance
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
        self.resolution_deg = res_deg
        self.tolerance = tol

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

    def generate_correction_matrix(self, par: AntennaPar, file_name: str):
        """
        Generates the correction factor matrix and saves it in a file

        Parameters:
            par (AntennaPar): set of antenna parameters to which calculate the
                correction factor
            file_name (str): name of file to which save the correction matrix
        """
        # Create antenna object
        azi = 0 # Antenna azimuth: 0 degrees for simplicity
        ele = 0 # Antenna elevation: 0 degrees as well
        self.antenna = AntennaBeamformingImt(par,azi,ele)

        # For co-channel beamforming
        # Correction factor numpy array
        correction_factor_co = np.zeros((len(self.phi_vals_deg),len(self.theta_vals_deg)))
        error_co = np.empty((len(self.phi_vals_deg),len(self.theta_vals_deg)), dtype=tuple)

        # Loop throug all the possible beams
        for phi_idx, phi in enumerate(self.phi_vals_deg):
            print('\n' + str(100*phi_idx/len(self.phi_vals_deg)) + '%')
            for theta_idx, theta in enumerate(self.theta_vals_deg):
                s = '\tphi = ' + str(phi) + ', theta = ' + str(theta)
                print(s)
                stdout.flush()
                correction_factor_co[phi_idx,theta_idx], error_co[phi_idx,theta_idx] = self.calculate_correction_factor(phi,theta,True)


        correction_factor_adj, error_adj = self.calculate_correction_factor(0,0,False)

        # Save in file
        self._save_files(correction_factor_co,
                         error_co,
                         correction_factor_adj,
                         error_adj,
                         par,
                         file_name)

    def calculate_correction_factor(self, phi_beam: float, theta_beam: float, c_chan: bool):
        """
        Calculates single correction factor

        Parameters:
            phi_beam (float): azimuth angle of beam pointing direction [deg]
            theta_beam (float): elevation angle of beam pointing direction [deg]
            c_chan (bool): if True, correction factor for whole antenna array
                is calculated. Otherwise, correction factor for a single
                element is calcualted

        Returns:
            correction_factor (float): correction factor value [dB]
            error (tuple): upper and lower error bounds [dB]
        """
        if c_chan:
            self.antenna.add_beam(phi_beam,theta_beam)
            beam = int(len(self.antenna.beams_list) - 1)
            int_f = lambda t,p: \
            np.power(10,self.antenna._beam_gain(np.rad2deg(p),np.rad2deg(t),beam)/10)*np.sin(t)
        else:
            int_f = lambda t,p: \
            np.power(10,self.antenna.element.element_pattern(np.rad2deg(p),np.rad2deg(t))/10)*np.sin(t)

        integral_val, err = dblquad(int_f,self.phi_min_rad,self.phi_max_rad,
                          lambda p: self.theta_min_rad,
                          lambda p: self.theta_max_rad,
                          epsabs=self.tolerance,
                          epsrel=0.0)

        correction_factor = -10*np.log10(integral_val/(4*np.pi))

        hig_bound = -10*np.log10((integral_val - err)/(4*np.pi))
        low_bound = -10*np.log10((integral_val + err)/(4*np.pi))

        return correction_factor, (low_bound,hig_bound)

    def _save_files(self, cf_co, err_co, cf_adj, err_adj, par, file_name):
        """
        Saves input correction factor and error values to npz file.
        Data is saved in an .npz file in a dict like data structure with the
        following keys:
            resolution (float): antenna array correction factor matrix angle
                resolution [deg]
            phi_range (tuple): range of beam pointing azimuth angle values
                [deg]
            theta_range (tuple): range of beam pointing elevation angle values
                [deg]
            correction_factor_co_channel (2D np.array): correction factor [dB]
                for the co-channel scenario (antenna array) for each of the phi
                theta pairs in phi_range and theta_range. Phi is associated
                with the lines and Theta is associated with the columns of the
                array
            error_co_channel (2D np.array of tuples): lower and upper bounds of
                calculated correction factors [dB], considering integration
                error
            correction_factor_adj_channel (float):correction factor for single
                antenna element
            error_adj_channel (tuple): lower and upper bounds [dB] of single
                antenna element correction factor
            parameters (AntennaPar): antenna parameters used in the
                normalization

        Parameters:
            cf_co (2D np.array): co-channel correction factor [dB]
            err_co (2D np.array): co-channel correction factor lower and upper
                bounds considering integration errors [dB]
            cf_adj (float): adjacent channel correction factor [dB]
            err_adj (float): adjancet channel correction factor lower and upper
                bounds considering integration errors [dB]
            par (AtennaPar): antenna parameters used in normalization
            file_name (str): name of file to which save normalization data
        """
        np.savez(file_name,
                 resolution = self.resolution_deg,
                 phi_range = (self.phi_min_deg, self.phi_max_deg),
                 theta_range = (self.theta_min_deg, self.theta_max_deg),
                 correction_factor_co_channel = cf_co,
                 error_co_channel = err_co,
                 correction_factor_adj_channel = cf_adj,
                 error_adj_channel = err_adj,
                 parameters = par)

if __name__ == '__main__':
    """
    Plots correction factor for horizontal and vertical planes.
    """
    import matplotlib.pyplot as plt
    import os

    # Create normalizer object
    resolution = 5
    tolerance = 1e-1
    norm = BeamformingNormalizer(resolution,tolerance)

    # Antenna parameters
    normalization = False
    norm_file = None
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
    par = AntennaPar(normalization,
                     norm_file,
                     element_pattern,
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
    file_name = 'main_test.npz'
    norm.generate_correction_matrix(par,file_name)
    data = np.load(file_name)
    correction_factor = data['correction_factor_co_channel']
    err_low, err_high = zip(*np.ravel(data['error_co_channel']))

    plt.plot(norm.phi_vals_deg,correction_factor,
             norm.phi_vals_deg,err_low,'r--',
             norm.phi_vals_deg,err_high,'r--')
    plt.xlim(-180,180)
    plt.ylabel(r"Correction factor [dB]")
    plt.xlabel(r"Azimuth angle $\phi$ [deg]")
    plt.title(r"Elevation angle $\theta$ = 90 deg")
    plt.show()
    data.close()

    # Set range of values & calculate correction factor
    norm = BeamformingNormalizer(resolution,tolerance)
    norm.phi_vals_deg = np.array([0])
    norm.generate_correction_matrix(par,file_name)
    data = np.load(file_name)
    correction_factor = data['correction_factor_co_channel']
    err_low, err_high = zip(*np.ravel(data['error_co_channel']))

    plt.plot(norm.theta_vals_deg,np.transpose(correction_factor),
             norm.theta_vals_deg,np.transpose(err_low),'r--',
             norm.theta_vals_deg,np.transpose(err_high),'r--')
    plt.xlim(0,180)
    plt.ylabel(r"Correction factor [dB]")
    plt.xlabel(r"Elevation angle $\theta$ [deg]")
    plt.title(r"Azimuth angle $\phi$ = 0 deg")
    plt.show()
    data.close()
    os.remove(file_name)

