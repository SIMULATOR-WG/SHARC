# -*- coding: utf-8 -*-
"""
Created on Fri Apr 14 14:13:58 2017

@author: Calil
"""

from sharc.antenna import Antenna

class AntennaImt(Antenna):
    """
    Implements a sectorized antenna, which corresponds to a single element of 
    an antenna array.
    
    Attributes
    ----------
        g_max (float): maximum gain of element
        theta_3db (float): vertical 3dB bandwidth of single element [degrees]
        phy_3db (float): horizontal 3dB bandwidth of single element [degrees]
        am (float): front-to-back ratio
        sla_v (float): element vertical sidelobe attenuation
    """
    
    def __init__(self,g_max: float, theta_3db: float, phy_3db: float, am: float, sla_v: float):
        """
        Constructs an AntennaImt object.
        
        Parameters
        ---------
            g_max (float): maximum gain of element
            theta_3db (float): vertical 3dB bandwidth of single element [degrees]
            phy_3db (float): horizontal 3dB bandwidth of single element [degrees]
            am (float): horizontal front-to-back ratio
            sla_v (float): vertical front-to-back ratio
        """
        self.__g_max = g_max
        self.__theta_3db = theta_3db
        self.__phy_3db = phy_3db
        self.__am = am
        self.__sla_v = sla_v
        
        # Create 
        super().__init__()
    
    @property
    def g_max(self):
        return self.__g_max
    
    @property
    def theta_3db(self):
        return self.__theta_3db
    
    @property
    def phy_3db(self):
        return self.__phy_3db
    
    @property
    def am(self):
        return self.__am
    
    @property
    def sla_v(self):
        return self.__sla_v
    
    def vertical_pattern(self,theta: float) -> float:
        """
        Calculates the vertical radiation pattern.
        
        Parameters
        ----------
            theta (float): elevation angle [degrees]
            
        Returns
        -------
            a_v (float): vertical radiation pattern gain value
        """
        pass
    
    def horizontal_pattern(self,phy: float) -> float:
        """
        Calculates the horizontal radiation pattern.
        
        Parameters
        ----------
            phy (float): azimuth angle [degrees]
            
        Returns
        -------
            a_h (float): horizontal radiation pattern gain value
        """
        pass
        
    def element_pattern(self, theta: float, phy: float) -> float:
        """
        Calculates the element radiation pattern gain.
        
        Parameters
        ----------
            theta (float): elevation angle [degrees]
            phy (float): azimuth angle [degrees]
            
        Returns
        -------
            a_e (float): element radiation pattern gain value
        """
        pass