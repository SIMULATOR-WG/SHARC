# -*- coding: utf-8 -*-
"""
Created on Mon May 15 12:51:48 2017

@author: LeticiaValle_Mac
"""

from sharc.propagation.propagation import Propagation 
from sharc.support.enumerations import StationType

import numpy as np
import scipy
import math
from scipy import special


class PropagationClutterLoss(Propagation):
    

    def __init__(self):
        super().__init__()
    
 
    def get_loss(self, *args, **kwargs) -> np.array:
        """
        Calculates clutter loss.
        
        Parameters
        ----------
            distance_3D (np.array) : 3D distances between stations
            distance_2D (np.array) : 2D distances between stations
            frequency (np.array) : center frequency [MHz]
            elevation (np.array) : elevation angles [deg]
            loc_percentage (np.array) : Percentage locations range [0, 1[
                                        "RANDOM" for random percentage
            station_type (StationType) : if type is IMT or FSS_ES, assume terrestrial 
                terminal within the clutter (ref § 3.2); otherwise, assume that
                one terminal is within the clutter and the other is a satellite, 
                aeroplane or other platform above the surface of the Earth.

        Returns
        -------
            array with clutter loss values with dimensions of frequency
        
        """
        f = kwargs["frequency"]
        theta = kwargs["elevation"]
        loc_per = kwargs["loc_percentage"]
        type = kwargs["station_type"]

        if loc_per == "RANDOM":
            p = np.random.random(f.shape)
        else:
            p = loc_per

        if type is StationType.IMT_BS or type is StationType.IMT_UE or type is StationType.FSS_ES:
            #Clutter Loss item 3.2 
            if "distance_2D" in kwargs:
                d = kwargs["distance_2D"]
            else:
                d = kwargs["distance_3D"]            
            
            Lt = 23.5 + 9.6*np.log10(f*1e-3)
            Ls = 32.98 + 23.9*np.log10(d*1e3) + 3*np.log10(f*1e-3)
            
            Q = np.sqrt(2)*scipy.special.erfcinv(2*(p))

            Lctt = -5*np.log10(10**(-0.2*Lt)+ 10**(-0.2*Ls)) - 6*Q 
 
        else:
            #Clutter Loss item 3.3 
            k1 = 93*(f*1e-3)**0.175
            A1 = 0.05
            
            y =np.sin(A1*(1 - (theta/90))+ math.pi*(theta/180))
            y1=np.cos(A1*(1 - (theta/90))+ math.pi*(theta/180))
            
            cot = (y1/y)                  
            Q = np.sqrt(2)*scipy.special.erfcinv(2*p)
            Lctt = (-k1*(np.log(1 - p))*cot)**(0.5*(90 - theta)/90) - 1 - 0.6*Q

        return Lctt


if __name__ == '__main__':
    import matplotlib.pyplot as plt

    #theta = np.array([90, 80, 70, 60, 50, 40, 30, 20, 15, 10, 5, 0])
    theta = np.array([90, 45, 30, 20 ])
    p = np.linspace(0, 1, 1001)
    freq = 27250*np.ones(theta.shape)
    
    cl = PropagationClutterLoss()
    clutter_loss = np.empty([len(theta), len(p)])
    
    for i in range(len(p)):
        clutter_loss[:,i] = cl.get_loss(frequency=freq,
                                        elevation=theta,
                                        loc_percentage=p[i],
                                        station_type=StationType.FSS_SS)
    
    fig = plt.figure(figsize=(8,6), facecolor='w', edgecolor='k')
    ax = fig.gca()

    for j in range(len(theta)):
        ax.plot(clutter_loss[j,:], 100*p, label="%i deg" % theta[j])
    
    plt.title("Cumulative distribution of clutter loss not exceeded for 27 GHz")
    plt.xlabel("clutter loss [dB]")
    plt.ylabel("percent of locations [%]")
    plt.xlim((-5, 20))
    plt.ylim((0, 100))                
    plt.legend(loc="lower right")
    plt.tight_layout()    
    plt.grid()

    plt.show()
        