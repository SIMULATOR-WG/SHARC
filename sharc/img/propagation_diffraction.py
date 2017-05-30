# -*- coding: utf-8 -*-


"""
Created on Tue Wen 17 10:15:31 2017

@author: LeticiaValle_Mac
"""

from sharc.propagation.propagation import Propagation 

import numpy as np
import math
 
class PropagationDiffraction(Propagation):
   
    """
    Implements the Diffraction model
    """
    def __init__(self):
        super(PropagationDiffraction, self).__init__()
        np.random.seed(0)

        #self.param = param
        #self.paramProp = paramProp
        self.propagation = PropagationDiffraction()
        
    def func_Gt(K,Bdft,Yt):

    
        Bt = Bdft*Yt
        
        if Bt>2:
            res_G = 17.6*(Bt - 1.1)**0.5 - 5*np.log10(Bt - 1.1) - 8
        else:
            res_G = 20*np.log(Bt + 0.1*Bt**3)
        
        K= 3
        
        if (res_G < (2 + 20*np.log10(K))):
            res_G = 2 + 20*np.log10(K)
        
        return res_G

    def func_Gr(K,Bdft,Yr):

        Br = Bdft*Yr
        
        if Br>2:
            res_G = 17.6*(Br - 1.1)**0.5 - 5*np.log10(Br - 1.1) - 8
        else:
            res_G = 20*np.log(Br + 0.1*Br**3)
        
       
        
        if (res_G < (2 + 20*np.log10(K))):
            res_G = 2 + 20*np.log10(K)
        
        return res_G


    def Jfunction(v):

        if (v<-0.78):
           v_result=0;
        
        else:
           v_result = 6.9 + 20*np.log10(np.sqrt((v - 0.1)**2 + 1) + v - 0.1)
        
        return v_result
    
    def Ld_difraction(d,f,a): # -*- coding: utf-8 -*-

        deltaN = 60
        hrs = 244
        hts = 280
        hte = 50
        hre = 50
        di= 2 
        hi = 6
        kbeta = 3
        polarization = 0 #horizontal
        er = 22.0
        sig = 0.003
        omega = 0
            
        k50 = 157/(157 - deltaN)
        ae = 6371*k50
        abeta = 6371*kbeta
        
        Ce = 1/ae
        lamb =(3*(10**8))/f
        
        
        #Implementing Ldsph
        if (a=='p'): 
           ap = ae
        if (a=='b'):
           ap = abeta
           
        m = (250*d**2)/(ap*(hte + hre))
        c = (hte - hre) / (hre + hre)
        b = 2*np.sqrt(((m+1)/3*m))*np.cos(math.pi + (1/3)*np.arccos((3*c/2)*np.sqrt((3*m)/(m+1)**3)))       
        dse1 = (d/2*(1+b))
        dse2 = (d - dse1)
        
        hse = ((hte - (500*dse1**2)/ap)*dse2 + (hre - (500*dse2**2)/ap)*dse1)/d
        hreq = 17.456*np.sqrt((dse1*dse2*lamb)/d)
        
        aem = 500*(d/(np.sqrt(hte)+np.sqrt(hre)))**2
        dlos = (np.sqrt(2*ap)*(np.sqrt(0.001*hte) + np.sqrt(0.001*hre)))
             
        if (d>dlos):
            adft = ap
        else:
             adft= aem
        
      
        for i in range(0,2,1): 
            if i == 0:
               er = 22.0
               sig = 0.003 
            
            if i == 1:
               er = 80.0
               sig = 5
               
            Kh = 0.036*(adft)**(-1/3)*((er - 1)**2 + (18*sig/f)**2)**(-1/4)
            Kv = Kh*(er**2 + (18*sig/f)**2)**0.5
        
            if (polarization == 0):  #horizontal polarization
               K = Kh
            else:
               K = Kv
        
            Bdft = (1 + 1.6*K**2 + 0.67*K**4)/(1 + 4.5*K**2 + 1.53*K**4)
            Yt = 0.9575*Bdft*(f/(adft**2))**(1/3)*hte
            Yr = 0.9575*Bdft*(f/(adft**2))**(1/3)*hre
            X = 21.88*Bdft*(f/(adft**2))**(1/3)*d
               
            if (X>=1.6):            
               Fx = 11 + 10*np.log10(X) - 17.6*X               
            else:
               Fx = -20*np.log10(X) - 5.6488*X**1.425
            
            if i == 0:                        
               Ldftland = -Fx - self.propagation.func_Gt(K,Bdft,Yt) - self.propagation.func_Gr(K,Bdft,Yr)
            if i == 1:                        
               Ldftsea = -Fx - self.propagation.func_Gt(K,Bdft,Yt) - self.propagation.func_Gr(K,Bdft,Yr)
             
            
        Ldft  = omega*Ldftsea + (1- omega)*Ldftland    
                                
        if (d>dlos):
            Ldsph = Ldft
       
        if(hse>hreq) or (Ldft<0):
            Ldsph = 0
        else:
            Ldsph = (1 - hse/hreq)*Ldft
        
        
        #Stim = max((hi + 500*Ce*di*(d - di) - hts)/di)
        #Srim = max((hi + 500*Ce*di(d - di) - hrs)/d - di)
        Stim = (hi + 500*Ce*di*(d - di) - hts)/di
        Srim = ((hi + 500*Ce*di*(d - di) - hrs)/(d - di))
    
        
        Dbp = (hrs - hts + Srim*d)/Stim +Srim
             
             
        if ((0.002*d)/(lamb*Dbp*(d - Dbp)<0)):      
           vb = (hts + Stim*Dbp - ((hts*(d - Dbp) + hrs*Dbp)/d))*(np.sqrt(0))
        else:
           vb = (hts + Stim*Dbp - ((hts*(d - Dbp) + hrs*Dbp)/d))*(np.sqrt((0.002*d)/(lamb*Dbp*(d - Dbp))))
           
            
        Luc = self.propagation.Jfunction(vb)  
       
        #Lbull = 2
        Lbull = Luc + (1 - np.exp(- Luc/6))*(10 + 0.02*d)
        Lbulla = Lbull
        Lbulls = Lbull
         
        Ld = Lbulla + max(Ldsph - Lbulls, 0)
        return Ld
    
    def get_loss(self, *args, **kwargs) -> np.array:
        
        d = np.asarray(kwargs["distance"])*(1e-3)  #Km
        f = np.asarray(kwargs["frequency"])*(1e-3)  #GHz
        a = 'p'
        f = 27
        
        if(a == 'p'):
            LdRes = self.propagation.Ld_difraction(d,f,a)    #Ld50
        
        if(a == 'b'):
            LdRes = self.propagation.Ld_difraction(d,f,a)  #Ldbeta
      
        
        return LdRes