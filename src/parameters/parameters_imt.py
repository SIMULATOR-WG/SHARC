# -*- coding: utf-8 -*-
"""
Created on Wed Feb 15 16:05:58 2017

@author: edgar
"""

class ParametersImt(object):

    ###########################################################################
    # Number of base stations per cluster
    num_base_stations = 19
    
    ###########################################################################
    # Number of clusters (should be 1 or 7)
    num_clusters = 1
    
    ###########################################################################
    # Configures static or dynamic positions for base stations
    static_base_stations = True
    
    ###########################################################################
    # Inter-site distance in macrocell network topology
    intersite_distance = 1500
    
    ###########################################################################
    # Defines if IMT service is the interferer or interfered-with service
    interfered_with = False
    
    ###########################################################################
    # IMT center frequency [MHz]
    frequency = 27250

    ###########################################################################
    # IMT bandwidth [MHz]
    bandwidth = 100    
    
    ###########################################################################
    # The load probability (or activity factor) models the statistical 
    # variation of the network load by defining the number of fully loaded
    # base stations that are simultaneously transmitting
    bs_load_probability = 0.5
    
    ###########################################################################
    # Base station transmit power [dBm]
    bs_tx_power = 37
    
    ###########################################################################
    # Base station height [m]
    bs_height = 30

    ###########################################################################
    # Number of UE that is allocated to each cell within to handover margin.
    # Remenber that in macrocell network each base station has 3 cells (sectors)
    ue_k = 30
    
    ###########################################################################
    # Multiplication factor that is used to ensure that the sufficient number
    # of UE's will distributed throughout ths system area such that the number
    # of K users is allocated to each cell. Normally, this values varies 
    # between 2 and 10 according to the user drop method
    ue_k_m = 10
    
    ###########################################################################
    # UE maximum transmit power [dBm]
    ue_tx_power = 37
    
    ###########################################################################
    # UE height [m]
    ue_height = 30