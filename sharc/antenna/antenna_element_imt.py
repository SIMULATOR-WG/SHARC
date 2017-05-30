# -*- coding: utf-8 -*-
"""
Created on Fri Apr 14 14:13:58 2017

@author: Calil
"""

import numpy as np

from sharc.support.named_tuples import AntennaPar

class AntennaElementImt(object):
    """
    Implements a single element of an IMT antenna array.
    
    Attributes
    ----------
        g_max (float): maximum gain of element
        theta_3db (float): vertical 3dB beamwidth of single element [degrees]
        phi_3db (float): horizontal 3dB beamwidth of single element [degrees]
        am (float): front-to-back ratio
        sla_v (float): element vertical sidelobe attenuation
    """
    
    def __init__(self,par: AntennaPar):
        """
        Constructs an AntennaElementImt object.
        
        Parameters
        ---------
            param (ParametersAntennaImt): antenna IMT parameters
        """        
        self.param = par
    
        self.__g_max = par.element_max_g
        self.__phi_3db = par.element_phi_3db
        self.__theta_3db = par.element_theta_3db
        self.__am = par.element_am
        self.__sla_v = par.element_sla_v
    
    @property
    def station_type(self):
        return self.__station_type
    
    @property
    def tx_or_rx(self):
        return self.__tx_or_rx
    
    @property
    def g_max(self):
        return self.__g_max
     
    @property
    def phi_3db(self):
        return self.__phi_3db
    
    @property
    def theta_3db(self):
        return self.__theta_3db
    
    @property
    def am(self):
        return self.__am
    
    @property
    def sla_v(self):
        return self.__sla_v
    
    def horizontal_pattern(self,phi: np.array) -> np.array:
        """
        Calculates the horizontal radiation pattern.
        
        Parameters
        ----------
            phi (np.array): azimuth angle [degrees]
            
        Returns
        -------
            a_h (np.array): horizontal radiation pattern gain value
        """
        return -1.0*np.minimum(12*(phi/self.phi_3db)**2,self.am)
    
    def vertical_pattern(self,theta: np.array) -> np.array:
        """
        Calculates the vertical radiation pattern.
        
        Parameters
        ----------
            theta (np.array): elevation angle [degrees]
            
        Returns
        -------
            a_v (np.array): vertical radiation pattern gain value
        """
        return -1.0*np.minimum(12*((theta-90.0)/self.theta_3db)**2,self.sla_v)
        
    def element_pattern(self, phi: np.array, theta: np.array) -> np.array:
        """
        Calculates the element radiation pattern gain.
        
        Parameters
        ----------
            theta (np.array): elevation angle [degrees]
            phi (np.array): azimuth angle [degrees]
            
        Returns
        -------
            gain (np.array): element radiation pattern gain value
        """
        att = -1.0*(self.horizontal_pattern(phi) + \
                    self.vertical_pattern(theta))
        gain = self.g_max - np.minimum(att,self.am)
        
        return gain