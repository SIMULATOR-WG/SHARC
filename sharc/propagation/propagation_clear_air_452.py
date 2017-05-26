# -*- coding: utf-8 -*-





#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 22 15:10:11 2017

@author: LeticiaValle_Mac
"""
from sharc.propagation.propagation import Propagation 
from sharc.propagation.propagation_gases_attenuation import PropagationGasesAttenuation
from sharc.propagation.propagation_duting_reflection import PropagationDutingReflection
from sharc.propagation.propagation_troposcatter import PropagationTropScatter

import numpy as np
 
class PropagationClearAir(Propagation):
    """
    Basic transmission loss due to free-space propagation and attenuation by atmospheric gases
    """
    def __init__(self):
        super(PropagationClearAir, self).__init__()
        np.random.seed(0)

        self.propagationAg = PropagationGasesAttenuation()
        self.propagationDuting = PropagationDutingReflection()
        self.propagationTropoScatter = PropagationTropScatter()
   
    def get_loss_Ag(self, *args, **kwargs) -> np.array:
    
        d = np.asarray(kwargs["distance"])*(1e-3)   #Km
        f = np.asarray(kwargs["frequency"])*(1e-3)  #GHz
        Ph = np.asarray(kwargs["atmospheric_pressure"])
        T = np.asarray(kwargs["air_temperature"])
        ro = np.asarray(kwargs["water_vapour"])
        Dlt = np.asarray(kwargs["Dlt"])
        Dlr = np.asarray(kwargs["Dlr"])
        Dct = np.asarray(kwargs["Dct"])
        Dcr = np.asarray(kwargs["Dcr"])
        
        Hts = np.asarray(kwargs["Hts"])
        Hrs = np.asarray(kwargs["Hrs"])
        Hte = np.asarray(kwargs["Hte"])
        Hre = np.asarray(kwargs["Hre"])
        thetaT = np.asarray(kwargs["theta_tx"])
        thetaR = np.asarray(kwargs["theta_rx"])
        N0 = np.asarray(kwargs["N0"])
        deltaN = np.asarray(kwargs["delta_N"])
        p = np.asarray(kwargs["percentage_p"])
       
        omega = np.asarray(kwargs["omega"])
        phi = np.asarray(kwargs["phi"])
        dtm = np.asarray(kwargs["dtm"])
        dlm = np.asarray(kwargs["dlm"])                                
        epsilon = np.asarray(kwargs["epsilon"])                                
        hm = np.asarray(kwargs["hm"])
        Gt = np.asarray(kwargs["tx_gain"])
        Gr = np.asarray(kwargs["rx_gain"])
        thetaT = np.asarray(kwargs["theta_tx"])
        thetaR = np.asarray(kwargs["theta_rx"])
        
        thetaJ = np.asarray(kwargs["thetaJ"])
        ep = np.asarray(kwargs["par_ep"])
        dsw = np.asarray(kwargs["dsw"])
        k = np.asarray(kwargs["k"])
        di = np.asarray(kwargs["dist_di"])
        hi = np.asarray(kwargs["hight_hi"])
        n = np.asarray(kwargs["eta"])
        Aht = np.asarray(kwargs["Aht"])
        Ahr = np.asarray(kwargs["Ahr"])
        
        C0 = np.asarray(kwargs["C0"])
        C1 = np.asarray(kwargs["C1"])
        C2 = np.asarray(kwargs["C2"])
        D1 = np.asarray(kwargs["D1"])
        D2 = np.asarray(kwargs["D2"])
        D3 = np.asarray(kwargs["D3"])
        
        #To find beta value                                   
        tau = 1 - np.exp(-(4.12*(10**-4)*dlm**2.41))
        mu1 = (10**(-dtm/(16 - 6.6*tau)) + (10**-(0.496 + 0.354*tau))**5)**0.2
        
        if (abs(phi)<=70):
           mu4 = 10**(-0.935 + 0.0176*abs(phi)*np.log10(mu1))
        else:
           mu4 = 10**(0.3*np.log10(mu1))
                
        if (abs(phi)<=70):
           beta = (10**(-0.015*abs(phi)+ 1.67)*mu1*mu4)
        else:
           beta = (4.17*mu1*mu4)    
        
        
        k50 = 157/(157 - deltaN)
        ae = 6371*k50
        Ce = 1/ae
        
        Stim = max((hi + 500*Ce*di*(d - di) - Hts)/di)
        Str = (Hrs - Hts)/d
        
        #Approximation to the inverse cumulative normal distribution function  
        Tx1 = np.sqrt((-2*np.log(p/100)))
        Ex1 = (((C2*Tx1+C1)*Tx1) +C0)/(((D3*Tx1 + D2)*Tx1 + D1)*Tx1 + 1)
        
        Tx2 = np.sqrt((-2*np.log(beta/100)))
        Ex2 = (((C2*Tx2+C1)*Tx2) +C0)/(((D3*Tx2 + D2)*Tx2 + D1)*Tx2 + 1)
        
        I1 = Ex1 - Tx1
        I2 = Ex2 - Tx2
        Fi = I1/ I2
        
        
        Fj = 1.0 - 0.5(1.0 + np.tanh(3.0*ep*(Stim - Str)/thetaJ))
        Fk = 1.0 - 0.5(1.0 + np.tanh(3.0*k*(d - dsw)/dsw))
        
        Lbfsg = self.propagationAg.get_loss_Ag(distance=d, frequency=f,atmospheric_pressure=Ph, air_temperature=T, water_vapour=ro)
        Lba =  self.propagationDuting.get_loss(distance=d, frequency=f,atmospheric_pressure=Ph, air_temperature= T, water_vapour=ro, Dlt = Dlt, Dlr=Dlr, Dct=Dct, Dcr=Dcr, Hts=Hts, Hrs=Hrs, Hte=Hte, Hre=Hre, theta_tx = thetaT, theta_rx = thetaR, N0 = N0, delta_N = deltaN, percentage_p = p, omega=omega, phi=phi, dtm=dtm, dlm=dlm, epsilon=epsilon, hm=hm)    #Ducting/layer reflection
        Lbs = self.propagationTropoScatter.get_loss(distance=d, frequency=f,atmospheric_pressure=Ph, air_temperature= T, water_vapour=ro, tx_gain = Gt, rx_gain = Gr, theta_tx = thetaT, theta_rx = thetaR, N0 = N0, delta_N = deltaN, percentage_p = p)       #Tropospheric scatter
        
        Esp = 2.6*(1 - np.exp(-0.1(Dlt + Dlr)))*np.log10(p/50)
        Esbeta = 2.6*(1 - np.exp(-0.1(Dlt + Dlr)))*np.log10(beta/50)
        Lb0p = Lbfsg + Esp
        Lb0beta = Lbfsg + Esbeta
        
        Ld50 = 0        #Difraction Loss not implemented
        Ldbeta = 0      #Difraction Loss not implemented
        
        Lbd50 = Lbfsg + Ld50
        Ldp = Ld50 + Fi*(Ldbeta - Ld50)
        Lbd = Lb0p + Ldp
        
        if (p < beta):
            Lminb0p = Lb0p + (1 - omega)*Ldp                   
        if (p >= beta):
            Lminb0p = Lbd50 + (Lb0beta + (1 - omega)*Ldp - Lbd50)*Fi
     
    
        Lminbap = n*np.log(np.exp(Lba/n) + np.exp(Lb0p/n))
       
        if (Lminbap > Lbd):
            Lbda = Lbd
        if (Lminbap <= Lbd):
            Lbda = Lminbap + (Lbd - Lminbap)*Fk
            
        Lbam = Lbda + (Lminb0p - Lbda)*Fj    
    
        Lb = -5*np.log10(10**(-0.2*Lbs) + 10**(-0.2*Lbam)) + Aht + Ahr
    
        Loss = Lb - Gt - Gr
        
        return Loss