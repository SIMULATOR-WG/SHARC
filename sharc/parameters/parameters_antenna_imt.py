# -*- coding: utf-8 -*-
"""
Created on Sat Apr 15 16:29:36 2017

@author: Calil
"""

from sharc.support.named_tuples import AntennaPar

class ParametersAntennaImt(object):
    """
    Defines parameters for antenna array.
    """
    
    def __init__(self):
        pass

    
    ###########################################################################
    # Named tuples which contain antenna types
    
    def get_antenna_parameters(self,sta_type: str, txrx: str)-> AntennaPar:
        if sta_type == "BS":
            if txrx == "TX":
                tpl = AntennaPar(self.bs_tx_element_max_g,\
                                       self.bs_tx_element_phi_3db,\
                                       self.bs_tx_element_theta_3db,\
                                       self.bs_tx_element_am,\
                                       self.bs_tx_element_sla_v,\
                                       self.bs_tx_n_rows,\
                                       self.bs_tx_n_columns,\
                                       self.bs_tx_element_horiz_spacing,\
                                       self.bs_tx_element_vert_spacing)
            elif txrx == "RX":
                tpl = AntennaPar(self.bs_rx_element_max_g,\
                                       self.bs_rx_element_phi_3db,\
                                       self.bs_rx_element_theta_3db,\
                                       self.bs_rx_element_am,\
                                       self.bs_rx_element_sla_v,\
                                       self.bs_rx_n_rows,\
                                       self.bs_rx_n_columns,\
                                       self.bs_rx_element_horiz_spacing,\
                                       self.bs_rx_element_vert_spacing)
        elif sta_type == "UE":
            if txrx == "TX":
                tpl = AntennaPar(self.ue_tx_element_max_g,\
                                       self.ue_tx_element_phi_3db,\
                                       self.ue_tx_element_theta_3db,\
                                       self.ue_tx_element_am,\
                                       self.ue_tx_element_sla_v,\
                                       self.ue_tx_n_rows,\
                                       self.ue_tx_n_columns,\
                                       self.ue_tx_element_horiz_spacing,\
                                       self.ue_tx_element_vert_spacing)
            elif txrx == "RX":
                tpl = AntennaPar(self.ue_rx_element_max_g,\
                                       self.ue_rx_element_phi_3db,\
                                       self.ue_rx_element_theta_3db,\
                                       self.ue_rx_element_am,\
                                       self.ue_rx_element_sla_v,\
                                       self.ue_rx_n_rows,\
                                       self.ue_rx_n_columns,\
                                       self.ue_rx_element_horiz_spacing,\
                                       self.ue_rx_element_vert_spacing)
                
        return tpl