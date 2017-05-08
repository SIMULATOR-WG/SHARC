# -*- coding: utf-8 -*-
"""
Created on Sat Apr 15 16:29:36 2017

@author: Calil
"""

from collections import namedtuple

class ParametersAntennaImt(object):
    """
    Defines parameters for antenna array.
    """
    
    __instance = None
    
    def __new__(cls):
        """
        This is the Singleton Pattern to ensure that this class will have only
        one instance
        """
        if ParametersAntennaImt.__instance is None:
            ParametersAntennaImt.__instance = object.__new__(cls)
        return ParametersAntennaImt.__instance
    
    ###########################################################################
    # Base station maximum element gain [dBi]
    bs_element_max_g = 5
    
    ###########################################################################
    # Base station horizontal 3dB beamwidth of single element [degrees]
    bs_element_phi_3db = 80
    
    ###########################################################################
    # Base station vertical 3dB beamwidth of single element [degrees]
    bs_element_theta_3db = 65
    
    ###########################################################################
    # Base station front to back ratio of single element [dB]
    bs_element_am = 30
    
    ###########################################################################
    # Base station element vertical sidelobe attenuation [dB]
    bs_element_sla_v = 30
    
    ###########################################################################
    # Base station number of rows in array
    bs_n_rows = 16
    
    ###########################################################################
    # Base station number of columns in array
    bs_n_columns = 16
    
    ###########################################################################
    # Base station array horizontal element spacing (d/lambda)
    bs_element_horiz_spacing = 0.5
    
    ###########################################################################
    # Base station array vertical element spacing (d/lambda)
    bs_element_vert_spacing = 0.5
    
    ###########################################################################
    # UE maximum element gain [dBi]
    ue_element_max_g = 5
    
    ###########################################################################
    # UE horizontal 3dB beamwidth of single element [degrees]
    ue_element_phi_3db = 80
    
    ###########################################################################
    # UE vertical 3dB beamwidth of single element [degrees]
    ue_element_theta_3db = 65
    
    ###########################################################################
    # UE front to back ratio of single element [dB]
    ue_element_am = 30
    
    ###########################################################################
    # UE element vertical sidelobe attenuation [dB]
    ue_element_sla_v = 30
    
    ###########################################################################
    # UE number of rows in array
    ue_n_rows = 4
    
    ###########################################################################
    # UE number of columns in array
    ue_n_columns = 2
    
    ###########################################################################
    # UE array horizontal element spacing (d/lambda)
    ue_element_horiz_spacing = 0.5
    
    ###########################################################################
    # UE array vertical element spacing (d/lambda)
    ue_element_vert_spacing = 0.5
    
    ###########################################################################
    # Named tuples which contain antenna types
    AntennaType = namedtuple("AntennaType","element_max_g element_phi_3db element_theta_3db element_am element_sla_v n_rows n_columns element_horiz_spacing element_vert_spacing")
    
    def get_antenna_parameters(self,sta_type: str)-> AntennaType:
        if sta_type == "BS":
            bs = self.AntennaType(self.bs_element_max_g,\
                                  self.bs_element_phi_3db,\
                                  self.bs_element_theta_3db,\
                                  self.bs_element_am,\
                                  self.bs_element_sla_v,\
                                  self.bs_n_rows,\
                                  self.bs_n_columns,\
                                  self.bs_element_horiz_spacing,\
                                  self.bs_element_vert_spacing)
            return bs
        elif sta_type == "UE":
            ue = self.AntennaType(self.ue_element_max_g,\
                                  self.ue_element_phi_3db,\
                                  self.ue_element_theta_3db,\
                                  self.ue_element_am,\
                                  self.ue_element_sla_v,\
                                  self.ue_n_rows,\
                                  self.ue_n_columns,\
                                  self.ue_element_horiz_spacing,\
                                  self.ue_element_vert_spacing)
            return ue
    