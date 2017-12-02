# -*- coding: utf-8 -*-
"""
Created on Tue Apr 25 15:05:38 2017

@author: LeticiaValle_Mac
"""
from sharc.propagation.propagation import Propagation

import numpy as np
import os

class PropagationGasesAttenuation(Propagation):
    """
    Basic transmission loss due to free-space propagation and attenuation by atmospheric gases
    """

    def get_loss(self, *args, **kwargs) -> np.array:

        d = np.asarray(kwargs["distance"]) #Km
        f = np.asarray(kwargs["frequency"]) #GHz

        Ph = np.asarray(kwargs["atmospheric_pressure"])
        T = np.asarray(kwargs["air_temperature"])
        ro = np.asarray(kwargs["water_vapour"])

        P=Ph
        theta=300/T
        e=ro*T/216.7
        p=P-e

        fo,a1,a2,a3,a4,a5,a6 = np.genfromtxt(os.path.join("..", "sharc", "propagation","data_oxygen.txt"), skip_header=1, unpack=True)
        fwv,b1,b2,b3,b4,b5,b6 = np.genfromtxt(os.path.join("..", "sharc", "propagation","data_water-vapour.txt"), skip_header=1, unpack=True)

        #Line strength Si
        Sio = np.array(a1)*(1e-7)*p*(theta**3.0)*np.exp(np.array(a2)*(1-theta))
        Siw = np.array(b1)*(1e-1)*e*(theta**3.5)*np.exp(np.array(b2)*(1-theta))

        #Line shape factor
        deltao = (np.array(a5)+np.array(a6)*theta)*(1e-4)*(p+e)*theta**(0.8)
        deltaw=0;

        deltafo=np.array(a3)*(1e-4)*(p*theta**(0.8-np.array(a4))+1.1*e*theta)
        deltafw=np.array(b3)*(1e-4)*(p*theta**np.array(b4)+np.array(b5)*e*theta**np.array(b6))

        part1Fio=(deltafo-deltao*(fo-f))/((fo-f)**2+deltafo**2)
        part2Fio=(deltafo-deltao*(fo+f))/((fo+f)**2+deltafo**2)
        Fio=(f/fo)*(part1Fio+part2Fio)

        part1Fiw=(deltafw-deltaw*(fwv-f))/((fwv-f)**2+deltafw**2)
        part2Fiw=(deltafw-deltaw*(fwv+f))/((fwv+f)**2+deltafw**2)
        Fiw=(f/fwv)*(part1Fiw+part2Fiw);

        #Dry air continuum N2Df ( N''D(f) )
        dis = 5.6*(1e-4)*(p+e)*theta**0.8;
        N2Df=f*p*theta**2*((6.14*(1e-5))/(dis*(1+(f/dis)**2))+(1.4*(1e-12)*p*theta**1.5)/(1+ 1.9*(1e-5)*f**1.5))

        N2f = sum(Fio*Sio)+sum(Fiw*Siw)+ N2Df

        gases_att = (0.1820*f*N2f);

        loss = gases_att*d

        return loss

