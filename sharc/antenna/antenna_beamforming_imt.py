# -*- coding: utf-8 -*-
"""
Created on Sat Apr 15 15:35:51 2017

@author: Calil
"""

import numpy as np

from sharc.antenna.antenna_imt import AntennaImt
from sharc.parameters.parameters_antenna_imt import ParametersAntennaImt

class AntennaBeamformingImt(AntennaImt):
    """
    Implements a an antenna array
    
    Attributes
    ----------
        position (np.array): x, y and z coordinates of array
        orientation (np.array): phi and theta of antenna orientation
            based on earth-centered coordinate system
        beam_orientation (np.array): phi_scan and theta_tilt of 
            antenna's beam
        gain (float): calculated antenna gain in given direction
        g_max (float): maximum gain of element
        theta_3db (float): vertical 3dB bandwidth of single element [degrees]
        phi_3db (float): horizontal 3dB bandwidth of single element [degrees]
        am (float): front-to-back ratio
        sla_v (float): element vertical sidelobe attenuation
        n_rows (int): number of rows in array
        n_cols (int): number of columns in array
        dh (float): horizontal element spacing over wavelenght (d/lambda)
        dv (float): vertical element spacing over wavelenght (d/lambda)
    """
    
    def __init__(self,param: ParametersAntennaImt, station_type: str):
        """
        Constructs an AntennaBeamformingImt object.
        
        Parameters
        ---------
            param (ParametersAntennaImt): antenna IMT parameters
            station_type (srt): type of station. Possible values are "BS" and
                "UE"
        """
        super().__init__(param, station_type)
        
        self.__position = np.array([0, 0, 0])

        self.__orientation = np.array([0, 0])
        
        self.__beam_orientation = np.array([0, 0])
        
        if station_type == "BS":
            self.__n_rows =param.bs_n_rows
            self.__n_cols =param.bs_n_columns
            self.__dh =param.bs_element_horiz_spacing
            self.__dv = param.bs_element_vert_spacing
        elif station_type == "UE":
            self.__n_rows =param.ue_n_rows
            self.__n_cols =param.ue_n_columns
            self.__dh =param.ue_element_horiz_spacing
            self.__dv = param.ue_element_vert_spacing
        
    def set_electrical_tilt(self,direct: np.array):
        """
        Set antenna electrical tilt to given earth-centered direction.
        
        Parameters
        ----------
            direct (np.array): coordinates of pointing target in 
                earth-coordinate system
        """
        pass
    
    @property
    def position(self):
        return self.__position
    
    @position.setter
    def position(self, pos):
        self.__position = pos
    
    @property
    def orientation(self):
        return self.__orientation
    
    @orientation.setter
    def orientation(self, ori):
        self.__orientation = ori
    
    @property
    def beam_orientation(self):
        return self.__beam_orientation
    
    @property
    def n_rows(self):
        return self.__n_rows
    
    @property
    def n_cols(self):
        return self.__n_cols
    
    @property
    def dh(self):
        return self.__dh
    
    @property
    def dv(self):
        return self.__dv
    
    def super_position_vector(self,phi: float, theta: float) -> np.array:
        """
        Calculates super position vector.
        All angles are taken according to the antenna coordinate system.
        
        Parameters
        ----------
            theta (float): elevation angle [degrees]
            phi (float): azimuth angle [degrees]
            
        Returns
        -------
            v_vec (np.array): superposition vector
        """
        r_phi = np.deg2rad(phi)
        r_theta = np.deg2rad(theta)
        
        n = np.arange(self.n_rows) + 1
        m = np.arange(self.n_cols) + 1
        
        exp_arg = (n[:,np.newaxis] - 1)*self.dv*np.cos(r_theta) + \
                  (m - 1)*self.dh*np.sin(r_theta)*np.sin(r_phi)
        
        v_vec = np.exp(2*np.pi*1.0j*exp_arg)
        
        return v_vec
        
    def weight_vector(self, phi_scan: float, theta_tilt: float) -> np.array:
        """
        Calculates super position vector.
        All angles are taken according to the antenna coordinate system.
        
        Parameters
        ----------
            phi_scan (float): electrical horizontal steering [degrees]
            theta_tilt (float): electrical down-tilt steering [degrees]
            
        Returns
        -------
            w_vec (np.array): weighting vector
        """
        r_phi = np.deg2rad(phi_scan)
        r_theta = np.deg2rad(theta_tilt)
        
        n = np.arange(self.n_rows) + 1
        m = np.arange(self.n_cols) + 1
        
        exp_arg = (n[:,np.newaxis] - 1)*self.dv*np.sin(r_theta) - \
                  (m - 1)*self.dh*np.cos(r_theta)*np.sin(r_phi)
        
        w_vec = (1/np.sqrt(self.n_rows*self.n_cols))*\
                np.exp(2*np.pi*1.0j*exp_arg)
        
        return w_vec
    
    def array_gain(self, v_vec: np.array, w_vec: np.array) -> float:
        """
        Calculates the array gain. Does not consider element gain.
        
        Parameters
        ----------
            v_vec (np.array): superposition vector
            w_vec (np.array): weighting vector
            
        Returns
        -------
            array_g (float): array gain
        """
        array_g = 10*np.log10(abs(np.sum(np.multiply(v_vec,w_vec)))**2)
        return array_g
        
    
    def beam_gain(self,phi: float, theta: float, phi_scan: float, theta_tilt: float) -> np.array:
        """
        Calculates gain for a single beam in a given direction.
        All angles are taken according to the antenna coordinate system.
        
        Parameters
        ----------
            phi (float): azimuth angle [degrees]
            theta (float): elevation angle [degrees]
            phi_scan (float): electrical horizontal steering [degrees]
            theta_tilt (float): electrical down-tilt steering [degrees]
            
        Returns
        -------
            gain (float): beam gain [dBi]
        """
        element_g = self.element_pattern(phi,theta)
        
        v_vec = self.super_position_vector(phi,theta)
        w_vec = self.weight_vector(phi_scan,theta_tilt)
        
        array_g = self.array_gain(v_vec,w_vec)
        
        self.gain = element_g + array_g
        
        return self.gain
    
    def set_gain_to(self, point: np.array):
        """
        Calculates gain to given point in the earth-centered coordinate system.
        
        Parameters
        ----------
            point (np.array): point to which gain is calculated in 
                earth-centered coordinate system
        """
        pass
    
    def convert_coordinate(self, point: np.array) -> np.array:
        """
        Converts direction in earth-centered coordinate system to antenna 
        coordinate system.
        
        Parameters
        ----------
            point (np.array): point in earth-centered coordinate system
            
        Returns
        -------
            coord_angles (np.array): azimuth and elevation angles
                based on antenna's coordination system
        """
        mov_point = point - self.__position
        
        print("mov_point = ",mov_point)
        
        r_orient = np.deg2rad(self.__orientation)
        r_phi = r_orient[0]
        r_theta = r_orient[1]
        
        rot_mtx_phi = np.array([[np.cos(r_phi), -np.sin(r_phi),   0],
                                [np.sin(r_phi),  np.cos(r_theta), 0],
                                [0,              0,               1]])
    
        print("rot_mtx_phi = ",rot_mtx_phi)
    
        rot_mtx_theta = np.array([[np.cos(r_theta),  0,   np.sin(r_theta)],
                                  [0,                1,                 0],
                                  [-np.sin(r_theta), 0,   np.cos(r_theta)]])
        
        print("rot_mtx_theta = ",rot_mtx_theta)    
        
        rot_mtx = np.dot(rot_mtx_phi,rot_mtx_theta)
        
        print("rot_mtx = ",rot_mtx)
        
        local_point = np.dot(rot_mtx,mov_point)
        
        print("local_point = ",local_point)
        
        r = np.linalg.norm(local_point)
        
        phi = np.arctan(local_point[1]/local_point[0])
        theta = np.arccos(local_point[2]/r)
        
        return np.rad2deg(np.array([phi, theta]))
    