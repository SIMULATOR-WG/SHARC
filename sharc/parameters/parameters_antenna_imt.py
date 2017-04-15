# -*- coding: utf-8 -*-
"""
Created on Sat Apr 15 16:29:36 2017

@author: Calil
"""

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
    bs_element_phy_3db = 80
    
    ###########################################################################
    # Base station vertical 3dB beamwidth of single element [degrees]
    bs_element_theta_3db = 80
    
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
    ue_element_phy_3db = 80
    
    ###########################################################################
    # UE vertical 3dB beamwidth of single element [degrees]
    ue_element_theta_3db = 80
    
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
    