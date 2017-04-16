# -*- coding: utf-8 -*-
"""
Created on Sat Apr 15 15:35:51 2017

@author: Calil
"""

import numpy as np

from sharc.antenna_imt import AntennaImt
from sharc.parameters.parameters_antenna_imt import ParametersAntennaImt

class AntennaBeamformingImt(AntennaImt):
    """
    Implements a an antenna array
    
    Attributes
    ----------
        gain (float): calculated antenna gain in given direction
        g_max (float): maximum gain of element
        theta_3db (float): vertical 3dB bandwidth of single element [degrees]
        phy_3db (float): horizontal 3dB bandwidth of single element [degrees]
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
    
    def super_position_vector(self,phy: float, theta: float) -> np.array:
        """
        Calculates super position vector.
        
        Parameters
        ----------
            theta (float): elevation angle [degrees]
            phy (float): azimuth angle [degrees]
            
        Returns
        -------
            v_vec (np.array): superposition vector
        """
        r_phy = np.deg2rad(phy)
        r_theta = np.deg2rad(theta)
        
        n = np.arange(self.n_rows) + 1
        m = np.arange(self.n_cols) + 1
        
        exp_arg = (n[:,np.newaxis] - 1)*self.dv*np.cos(r_theta) + \
                  (m - 1)*self.dh*np.sin(r_theta)*np.sin(r_phy)
        
        v_vec = np.exp(2*np.pi*1.0j*exp_arg)
        
        return v_vec
        
    def weight_vector(self, phy_scan: float, theta_tilt: float) -> np.array:
        """
        Calculates super position vector.
        
        Parameters
        ----------
            phy_scan (float): electrical horizontal steering [degrees]
            theta_tilt (float): electrical down-tilt steering [degrees]
            
        Returns
        -------
            w_vec (np.array): weighting vector
        """
        r_phy = np.deg2rad(phy_scan)
        r_theta = np.deg2rad(theta_tilt)
        
        n = np.arange(self.n_rows) + 1
        m = np.arange(self.n_cols) + 1
        
        exp_arg = (n[:,np.newaxis] - 1)*self.dv*np.sin(r_theta) - \
                  (m - 1)*self.dh*np.cos(r_theta)*np.sin(r_phy)
        
        w_vec = (1/np.sqrt(self.n_rows*self.n_cols))*\
                np.exp(2*np.pi*1.0j*exp_arg)
        
        return w_vec
    
    def array_gain(self, v_vec: np.array, w_vec: np.array) -> float:
        """
        Calculates the array gain. Does not consider element gain,
        
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
        
    
    def beam_gain(self,phy: float, theta: float, phy_scan: float, theta_tilt: float) -> np.array:
        """
        Calculates gain for a single beam in a given direction.
        
        Parameters
        ----------
            phy (float): azimuth angle [degrees]
            theta (float): elevation angle [degrees]
            phy_scan (float): electrical horizontal steering [degrees]
            theta_tilt (float): electrical down-tilt steering [degrees]
            
        Returns
        -------
            beam_gain (float): beam gain [dB]
        """
        pass