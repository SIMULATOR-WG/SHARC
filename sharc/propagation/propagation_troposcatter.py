"""
Created on Tue Mai 08 12:05:38 2017

@author: LeticiaValle_Mac
"""
from sharc.propagation.propagation import Propagation
from sharc.propagation.propagation_gases_attenuation import PropagationGasesAttenuation

import numpy as np

class PropagationTropScatter(Propagation):
    """
    Basic transmission loss due to free-space propagation and attenuation by atmospheric gases
    """

    def __init__(self, random_number_gen: np.random.RandomState, propagation_ag: PropagationGasesAttenuation):
        super().__init__(random_number_gen)

        #self.param = param
        #self.paramProp = paramProp
        self.propagation = propagation_ag

    def get_loss(self, *args, **kwargs) -> np.array:
        #loss = self.propagation.get_loss(distance=d, frequency=self.param.frequency)

        d = np.asarray(kwargs["distance"])*(1e-3)   #Km
        f = np.asarray(kwargs["frequency"])*(1e-3)  #GHz
        Ph = np.asarray(kwargs["atmospheric_pressure"])
        T = np.asarray(kwargs["air_temperature"])
        ro = np.asarray(kwargs["water_vapour"])

        Gt = np.asarray(kwargs["tx_gain"])
        Gr = np.asarray(kwargs["rx_gain"])
        thetaT = np.asarray(kwargs["theta_tx"])
        thetaR = np.asarray(kwargs["theta_rx"])
        No = np.asarray(kwargs["N0"])
        deltaN = np.asarray(kwargs["delta_N"])
        p = np.asarray(kwargs["percentage_p"])
        number_of_sectors = kwargs["number_of_sectors"]

        loss_Ag = self.propagation.get_loss(distance=d, frequency=f,atmospheric_pressure=Ph, air_temperature=T, water_vapour=ro)

        Lf = 25*np.log10(f) - 2.5*(np.log10(f/2))**2      #Frequency dependent loss (dB)
        Lc = 0.051*np.exp(0.055*(Gt + Gr))                #Aperture to medium coupling loss (dB)

        #Definition of the angular distance (mrad)
        k50 = 157/(157 - deltaN)
        Ae = 6371*k50
        teta = d*(10**3)/Ae +thetaT + thetaR

        if number_of_sectors > 1:
            d = np.repeat(d, number_of_sectors, 1)
            teta = np.repeat(teta, number_of_sectors, 1)
            loss_Ag = np.repeat(loss_Ag, number_of_sectors, 1)

        loss = 190 + Lf + 20*np.log10(d) +0.573*teta - 0.15*No + Lc + loss_Ag - 10.1*(-np.log10(p/50))**0.7

        return loss
