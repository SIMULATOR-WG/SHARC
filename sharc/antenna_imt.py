# -*- coding: utf-8 -*-
"""
Created on Fri Apr 14 14:13:58 2017

@author: Calil
"""

from sharc.antenna import Antenna
from sharc.parameters.parameters_antenna_imt import ParametersAntennaImt

class AntennaImt(Antenna):
    """
    Implements a sector antenna, which corresponds to a single element of 
    an antenna array.
    
    Attributes
    ----------
        gain (float): calculated antenna gain in given direction
        g_max (float): maximum gain of element
        theta_3db (float): vertical 3dB beamwidth of single element [degrees]
        phy_3db (float): horizontal 3dB beamwidth of single element [degrees]
        am (float): front-to-back ratio
        sla_v (float): element vertical sidelobe attenuation
    """
    
    def __init__(self,param: ParametersAntennaImt, station_type: str):
        """
        Constructs an AntennaImt object.
        
        Parameters
        ---------
            param (ParametersAntennaImt): antenna IMT parameters
            station_type (srt): type of station. Possible values are "BS" and
                "UE"
        """
        self.param = param
        
        self.__station_type = station_type
        
        if station_type == "BS":
            self.__g_max = param.bs_element_max_g
            self.__phy_3db = param.bs_element_phy_3db
            self.__theta_3db = param.bs_element_theta_3db
            self.__am = param.bs_element_am
            self.__sla_v = param.bs_element_sla_v
        elif station_type == "UE":
            self.__g_max = param.ue_element_max_g
            self.__phy_3db = param.ue_element_phy_3db
            self.__theta_3db = param.ue_element_theta_3db
            self.__am = param.ue_element_am
            self.__sla_v = param.ue_element_sla_v
        
        # Call for super class constructor 
        super().__init__()
    
    @property
    def station_type(self):
        return self.__station_type
    
    @property
    def g_max(self):
        return self.__g_max
    
        
    @property
    def phy_3db(self):
        return self.__phy_3db
    
    @property
    def theta_3db(self):
        return self.__theta_3db
    
    @property
    def am(self):
        return self.__am
    
    @property
    def sla_v(self):
        return self.__sla_v
    
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
        return -1.0*min(12*(phy/self.phy_3db)**2,self.am)
    
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
        return -1.0*min(12*((theta-90.0)/self.theta_3db)**2,self.sla_v)
        
    def element_pattern(self, phy: float, theta: float) -> float:
        """
        Calculates the element radiation pattern gain.
        
        Parameters
        ----------
            theta (float): elevation angle [degrees]
            phy (float): azimuth angle [degrees]
            
        Returns
        -------
            gain (float): element radiation pattern gain value
        """
        att = -1.0*(self.horizontal_pattern(phy) + \
                    self.vertical_pattern(theta))
        self.gain = self.g_max - min(att,self.am)
        
        return self.gain