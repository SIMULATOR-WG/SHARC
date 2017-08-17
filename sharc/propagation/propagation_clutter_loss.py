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
            distance_3D (np.array) : 3D distances between stations [m]
            distance_2D (np.array) : 2D distances between stations [m]
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
        loc_per = kwargs["loc_percentage"]
        type = kwargs["station_type"]

        if loc_per == "RANDOM":
            p = np.random.random(f.shape)
        else:
            p = loc_per*np.ones(f.shape)

        if type is StationType.IMT_BS or type is StationType.IMT_UE or type is StationType.FSS_ES:
            # Clutter Loss item 3.2 
            # Statistical clutter loss model for terrestrial paths
            if "distance_2D" in kwargs:
                d = kwargs["distance_2D"]
            else:
                d = kwargs["distance_3D"]            
            
            # minimum path length for the correction to be applied at only one end of the path
            idx = np.where(d > 250)[0]
            
            if len(idx):
                Lt = 23.5 + 9.6*np.log10(f[idx]*1e-3)
                Ls = 32.98 + 23.9*np.log10(d[idx]*1e-3) + 3*np.log10(f[idx]*1e-3)
                
                Q = np.sqrt(2)*scipy.special.erfcinv(2*(p[idx]))
    
                loss = -5*np.log10(10**(-0.2*Lt)+ 10**(-0.2*Ls)) - 6*Q 
    
                Lctt = np.zeros(d.shape)
                Lctt[idx] = loss
            else:
                Lctt = np.zeros(d.shape)
 
        else:
            # Clutter Loss item 3.3 
            # Earth-space and Aeronautical statistical clutter loss
            theta = kwargs["elevation"]

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
    #from cycler import cycler

    ###########################################################################
    # Earth-space and Aeronautical statistical clutter loss
    theta = np.array([90, 80, 70, 60, 50, 40, 30, 20, 15, 10, 5, 0])
    #theta = np.array([90, 45, 30, 20 ])
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
    #ax.set_prop_cycle( cycler('color', ['k', 'r', 'b', 'g']) )

    for j in range(len(theta)):
        ax.plot(clutter_loss[j,:], 100*p, label="%i deg" % theta[j], linewidth=2)
    
    plt.title("Cumulative distribution of clutter loss not exceeded for 27 GHz")
    plt.xlabel("clutter loss [dB]")
    plt.ylabel("percent of locations [%]")
    plt.xlim((-10, 70))
    plt.ylim((0, 100))                
    plt.legend(loc="lower right")
    plt.tight_layout()    
    plt.grid()

    ###########################################################################
    # Statistical clutter loss model for terrestrial paths    
    
    distance = np.linspace(250, 100000, 100000)
    frequency = np.array([2, 3, 6, 16, 40, 67])*1e3
    
    clutter_loss_ter = np.empty([len(frequency), len(distance)])
    
    for i in range(len(frequency)):
            clutter_loss_ter[i,:] = cl.get_loss(frequency = frequency[i] * np.ones(distance.shape),
                                            distance_2D = distance,
                                            loc_percentage = 0.5,
                                            station_type = StationType.FSS_ES)           
    
    fig = plt.figure(figsize=(8,6), facecolor='w', edgecolor='k')
    ax = fig.gca()
    #ax.set_prop_cycle( cycler('color', ['k', 'r', 'b', 'g']) )

    for j in range(len(frequency)):
        freq = frequency[j]*1e-3
        ax.semilogx(distance*1e-3, clutter_loss_ter[j,:], label="%i GHz" % freq, linewidth=2)
    
    plt.title("Median clutter loss for terrestrial paths")
    plt.xlabel("Distance [km]")
    plt.ylabel("Median clutter loss [dB]")
    plt.xlim((0.2, 100))
    plt.ylim((15, 45))      
    ax.set_xticks(np.linspace(0.2, 1, 8).tolist() + np.linspace(2, 10, 9).tolist() + np.linspace(20, 100, 9).tolist())          
    plt.legend(loc="lower right")
    plt.tight_layout()    
    plt.grid()        
        
    plt.show()
        