# -*- coding: utf-8 -*-
"""
Created on Mon May 15 12:51:48 2017

@author: LeticiaValle_Mac
"""

from sharc.propagation.propagation import Propagation 
from sharc.parameters.parameters_imt import ParametersImt
from sharc.parameters.parameters_propagation import ParametersPropagation

import numpy as np
 
import scipy
import math
from scipy.stats import norm
from scipy.special import erf

class PropagationClutterLoss(Propagation):
 
  def get_loss(self, *args, **kwargs) -> np.array:
      
      f = np.asarray(kwargs["frequency"])
      p = np.asarray(kwargs["percentage_p"])
      d = np.asarray(kwargs["dist"])
      theta = np.asarray(kwargs["elevation_angle_facade"])
      P = np.asarray(kwargs["probability_loss_notExceeded"])
      r = np.asarray(kwargs["coeff_r"])
      s = np.asarray(kwargs["coeff_s"])
      t = np.asarray(kwargs["coeff_t"])
      u = np.asarray(kwargs["coeff_u"])
      v = np.asarray(kwargs["coeff_v"])
      w = np.asarray(kwargs["coeff_w"])
      x = np.asarray(kwargs["coeff_x"])
      y = np.asarray(kwargs["coeff_y"])
      z = np.asarray(kwargs["coeff_z"])
      
     #Building Entry Loss   
      Lh = r + s*np.log10(f) + t*(np.log10(f))**2
      Le = 0.212*abs(theta)  
    
      C = -3
      mu1 = Lh + Le
      mu2 = w + x*np.log10(f)
      sig1 = u + v*np.log10(f)
      sig2 = y + z*np.log10(f)
    
             
      Ap = scipy.stats.norm.ppf(P)*sig1 + mu1
      Bp = scipy.stats.norm.ppf(P)*sig2 + mu2
 
                             
      Lbel = 10*np.log10(10**(0.1*Ap) +10**(0.1*Bp) + (10**(0.1*C)))
      
      
      
      #Clutter Loss item 3.1 
      Lt = 23.5 + 9.6*np.log10(f)
      Ls = 32.98 + 23.9*np.log10(d) + 3*np.log10(f)

      Q = (scipy.stats.norm.ppf(p/100))**-1
    
      Lctt = -5*np.log10(10**(-0.2*Lt)+ 10**(-0.2*Ls)) - 6*Q
     
      return (Lbel + Lctt)     
 
  
      
                
                