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
    Implements the ABG loss model according to the article "Propagation Path 
    Loss Models for 5G Urban Microand Macro-Cellular Scenarios"
    """
    
    
    def __init__(self):
        super().__init__()
        self.alpha = 3.4
        self.beta = 19.2
        self.gamma = 2.3
        self.building_loss = 20
        self.shadowing_sigma_dB = 6.5
    
    
    def get_loss(self, *args, **kwargs) -> np.array:
        """
        Calculates path loss for LOS and NLOS cases with respective shadowing
        (if shadowing is to be added)
        
        Parameters
        ----------
            distance_2D (np.array) : distances between stations [m]
            frequency (np.array) : center frequencie [MHz]
            line_of_sight_prob (float) : probability of LOS
            alpha (float): captures how the PL increases as the distance increases
            beta (float): floating offset value in dB
            gamma(float): captures the PL variation over the frequency
            shadowing (bool) : standard deviation value

        Returns
        -------
            array with path loss values with dimensions of distance_2D
        
        """
        d = kwargs["distance_2D"]
        f = kwargs["frequency"]
        p_los = kwargs["line_of_sight_prob"]

        building_loss = np.random.choice([0, self.building_loss], d.shape, p=[p_los, 1-p_los])
        
        if "alpha" in kwargs:
            self.alpha = kwargs["alpha"]
            
        if "beta" in kwargs:
            self.beta = kwargs["beta"]
            
        if "gamma" in kwargs:
            self.gamma = kwargs["gamma"]
            
        if "shadowing" in kwargs:
            std = kwargs["shadowing"] 
        else:
            std = False
        
    
        if std:
            shadowing = np.random.normal(0, self.shadowing_sigma_dB, d.shape)
        else:
            shadowing = 0
            
       
        loss = 10*self.alpha*np.log10(d) + self.beta + 10*self.gamma*np.log10(f*1e-3) + \
                shadowing + building_loss
        
        return loss

   
  
if __name__ == '__main__':
        
    ###########################################################################
    # Print path loss for ABG and Free Space models
    from propagation_free_space import PropagationFreeSpace
    
    ABG = PropagationABG()
    
    shadowing_std = 0
    distance = np.linspace(1, 10000, num=10000)[:,np.newaxis]
    frequency = 27000*np.ones(distance.shape)
    alphaUMA = 3.4
    betaUMA = 19.2
    gammaUMA = 2.3
    x_sfUMA = 0
    
    alphaUMI = 3.53
    betaUMI = 22.4
    gammaUMI = 2.13
    x_sfUMI = 0
    
    lossUMA_NLOS = ABG.get_loss(distance_2D = distance, frequency = frequency, alpha = alphaUMA, beta = betaUMA, gamma = gammaUMA)
    lossUMI_NLOS = ABG.get_loss(distance_2D = distance, frequency = frequency, alpha = alphaUMI, beta = betaUMI, gamma = gammaUMI)
    loss_fs = PropagationFreeSpace().get_loss(distance_2D=distance, frequency=frequency)
    
    fig = plt.figure(figsize=(8,6), facecolor='w', edgecolor='k')
    ax = fig.gca()
    ax.set_prop_cycle( cycler('color', ['r', 'g', 'b', 'y']) )

    ax.semilogx(distance, lossUMA_NLOS, label="UMa NLOS scenario")
    ax.semilogx(distance, lossUMI_NLOS, label="Umi Street Canyon NLOS scenario")
    ax.semilogx(distance, loss_fs, label="Free space")
        
    plt.title("ABG path loss model")
    plt.xlabel("distance [m]")
    plt.ylabel("path loss [dB]")
    plt.xlim((0, distance[-1,0]))
    plt.ylim((60, 200))                
    plt.legend(loc="upper left")
    plt.tight_layout()    
    plt.grid()

    plt.show()
    