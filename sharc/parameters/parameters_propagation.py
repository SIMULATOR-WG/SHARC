# -*- coding: utf-8 -*-
"""
Created on Fri May  5 16:36:47 2017

@author: edgar
"""


class ParametersPropagation(object):

    __instance = None

    def __new__(cls):
        """
        This is the Singleton Pattern to ensure that this class will have only
        one instance
        """
        if ParametersPropagation.__instance is None:
            ParametersPropagation.__instance = object.__new__(cls)
        return ParametersPropagation.__instance


    ###########################################################################
    # Total air pressure in hPa
    atmospheric_pressure = 1013
    
    ###########################################################################
    # Temperature in Kelvin
    air_temperature = 288
    
    ###########################################################################
    # water vapour concentration (g/m^3)
    water_vapour = 3
    
    ###########################################################################
    # Transmit antena gain (dB)
    tx_gain = 10
    
    ###########################################################################
    # Receive antena gain (dB)
    rx_gain = 10
    
    ###########################################################################
    #Transmit horizon elevation (mrad)
    theta_tx = 20
    
    ###########################################################################
    #Receive horizon elevation (mrad)
    theta_rx = 20
   
    ###########################################################################
    #Sea-level surface refractivity (use the map)
    N0 = 355
    
    ###########################################################################
    #Average radio-refractive (use the map)
    delta_N = 60
    
    ###########################################################################
    #percentage p
    percentage_p = 40
    
    ###########################################################################
    #distance from the transmit antennas to their respective horizons (km
    Dlt = 30
    
    ###########################################################################
    #distance from the receive antennas to their respective horizons (km)
    Dlr = 10 
    
    ###########################################################################
    #Distance over land from the transmit and receive antennas to the coast (km)
    Dct = 10
    
    ###########################################################################
    #Distance over land from the transmit and receive antennas to the coast (km)
    Dcr = 10 
    
    ###########################################################################
    #Antenna centre height above mean sea level (m)
    Hts = 280
    
    ###########################################################################
    #Antenna centre height above mean sea level (m)
    Hrs = 244
    
    ###########################################################################
    ##Effective height of interfering antenna (m)
    Hte = 50   
    
    ###########################################################################
    #Effective height of interfered-with antenna (m)
    Hre = 50       
    
    ###########################################################################
    #Fraction of the total path over water, 0 for totally overland paths
    omega = 0 
    
    ###########################################################################
    #path centre latitude (degrees)
    phi = 60 
    
    ###########################################################################
    #longest continuous land (inland + coastal) section of the great-circle path (km)
    dtm = 0.8 
    
    ###########################################################################
    #longest continuous inland section of the great-circle path (km)  
    dlm = 0.8 
    
    ###########################################################################
    #e = 3.5
    epsilon = 3.5
    
    ###########################################################################
    #nao sei o que é, pag. 18 eq (56), ITU -452
    hm = 15
    
    ###########################################################################
    #Elevation angle of the path at the building façade (degrees)
    elevation_angle_facade = 0
    
    ###########################################################################
    #Probability that loss is not exceeded
    probability_loss_notExceeded = 0.9
    
    ###########################################################################
    #Model coeficients
    coeff_r = 12.64
    coeff_s = 3.72
    coeff_t = 0.96
    coeff_u = 9.6
    coeff_v = 2
    coeff_w = 9.1
    coeff_x = -3
    coeff_y = 4.5
    coeff_z = -2
    
    ###########################################################################
    #Distance in km
    dist = 10