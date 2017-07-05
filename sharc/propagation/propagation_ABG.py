# -*- coding: utf-8 -*-

# -*- coding: utf-8 -*-
"""
Created on Tue Jul  4 11:57:41 2017

@author: LeticiaValle_Mac
"""

from sharc.propagation.propagation import Propagation 

import numpy as np
import matplotlib.pyplot as plt
from cycler import cycler
 
class PropagationABG(Propagation):
    """
    Implements the ABG loss model according to the article "Propagation Path Loss Models for 5G Urban Microand
    Macro-Cellular Scenarios"
    """
    
    
    def __init__(self):
        super(Propagation, self).__init__()    
    
    
    def get_loss(self, *args, **kwargs) -> np.array:
        """
        Calculates path loss for LOS and NLOS cases with respective shadowing
        (if shadowing is to be added)
        
        Parameters
        ----------
            distance (np.array) : distances between stations (m)
            frequency (float) : center frequencie [MHz]
            shadowing (float) : SF term in dB
            alpha (float): captures how the PL increases as the distance increases
            beta (float): floating offset value in dB
            gama(float): captures the PL variation over the frequency
            shadowing (bool) : if shadowing should be added or not
        Returns
        -------
            array with path loss values with dimensions of distance_2D
        
        """
        d = kwargs["distance"]
        f = kwargs["frequency"]*(1e-3) #converts to GHz
        alpha = kwargs["ABG_alpha"]
        beta = kwargs["ABG_beta"]
        gamma = kwargs["ABG_gamma"]
        x_sf = kwargs["shadowing"] 
        
    
        if x_sf:
            shadowing_std = np.random.normal(0, x_sf)
        else:
            shadowing_std = 0
            
       
        loss = 10*alpha*np.log10(d) + beta + 10*gamma*np.log10(f) + shadowing_std
        
        return loss

   
  
if __name__ == '__main__':
    
   
        
    ###########################################################################
    # Print path loss for ABG and Free Space models
    from propagation_free_space import PropagationFreeSpace
    
    ABG = PropagationABG()
    
    shadowing_std = 0
    distance = np.linspace(1, 10000, num=10000)[:,np.newaxis]
    frequency = 27000
    alphaUMA = 3.4
    betaUMA = 19.2
    gammaUMA = 2.3
    x_sfUMA = 6.5
    
    alphaUMI = 3.53
    betaUMI = 22.4
    gammaUMI = 2.13
    x_sfUMI = 7.82
    
    lossUMA_NLOS = ABG.get_loss(distance = distance, frequency = frequency, ABG_alpha = alphaUMA, ABG_beta = betaUMA, ABG_gamma = gammaUMA, shadowing = x_sfUMA)
    lossUMI_NLOS = ABG.get_loss(distance = distance, frequency = frequency, ABG_alpha = alphaUMI, ABG_beta = betaUMI, ABG_gamma = gammaUMI, shadowing = x_sfUMI)
    loss_fs = PropagationFreeSpace().get_loss(distance=distance, frequency=frequency)
    
    fig = plt.figure(figsize=(8,6), facecolor='w', edgecolor='k')
    ax = fig.gca()
    ax.set_prop_cycle( cycler('color', ['r', 'g', 'b', 'y']) )

    ax.plot(distance, lossUMA_NLOS, label="UMa NLOS scenario")
    ax.plot(distance, lossUMI_NLOS, label="Umi Street Canyon NLOS scenario")
    ax.plot(distance, loss_fs, label="Free space")
        
    plt.title("ABG path loss model")
    plt.xlabel("distance [m]")
    plt.ylabel("path loss [dB]")
    plt.xlim((0, distance[-1,0]))
    plt.ylim((60, 200))                
    plt.legend(loc="upper left")
    plt.tight_layout()    
    plt.grid()

    plt.show()
    