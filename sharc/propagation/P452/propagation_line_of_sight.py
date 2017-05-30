# -*- coding: utf-8 -*-

"""
Created on Tue Wen 17 10:15:31 2017

@author: LeticiaValle_Mac
"""


from sharc.propagation.P452.propagation_gases_attenuation import PropagationGasesAttenuation
from sharc.propagation.propagation import Propagation 

import numpy as np
 
class PropagationLineOfSight(Propagation):
    """
    Implements the Free Space propagation model
    """
    
    def __init__(self):
        super(PropagationLineOfSight, self).__init__()
        np.random.seed(0)

        #self.param = param
        #self.paramProp = paramProp
        self.propagation = PropagationGasesAttenuation()
        
    def get_loss(self, *args, **kwargs) -> np.array:
        d = np.asarray(kwargs["distance"])*(1e-3)  #Km
        f = np.asarray(kwargs["frequency"])*(1e-3)  #GHz
        Ph = np.asarray(kwargs["atmospheric_pressure"])
        T = np.asarray(kwargs["air_temperature"])
        ro = np.asarray(kwargs["water_vapour"])
        
        loss_Ag = self.propagation.get_loss_Ag(distance=d, frequency=f,atmospheric_pressure=Ph, air_temperature=T, water_vapour=ro)
        loss = 20*np.log10(d) + 20*np.log10(f) + 92.5  + loss_Ag
                          
        return loss