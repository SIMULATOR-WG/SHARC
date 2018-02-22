#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 24 12:10:11 2017

@author: LeticiaValle_Mac
"""

from sharc.propagation.propagation_gases_attenuation import PropagationGasesAttenuation
from sharc.propagation.propagation import Propagation

import math
import numpy as np


class PropagationDiffraction(Propagation):
    """
    Implements the Diffraction loss model
    """

    def __init__(self, random_number_gen, propagation_ag: PropagationGasesAttenuation):
        super().__init__(random_number_gen)

        self.propagationAg = propagation_ag

    def func_Gt(self, K,Bdft,Yt):

        Bt = Bdft*Yt

        if Bt>2:
            res_G = 17.6*(Bt - 1.1)**0.5 - 5*np.log10(Bt - 1.1) - 8
        else:
            res_G = 20*np.log10(Bt + 0.1*Bt**3)


        if (res_G < (2 + 20*np.log10(K))):
            res_G = 2 + 20*np.log10(K)

        return res_G

    def func_Gr(self,K,Bdft,Yr):

        Br = Bdft*Yr

        if Br>2:
            res_G = 17.6*(Br - 1.1)**0.5 - 5*np.log10(Br - 1.1) - 8
        else:
            res_G = 20*np.log10(Br + 0.1*Br**3)



        if (res_G < (2 + 20*np.log10(K))):
            res_G = 2 + 20*np.log10(K)

        return res_G

    def Jfunction(self,v):
        if (v<-0.78):
           v_result=0;

        else:
           v_result = 6.9 + 20*np.log10(np.sqrt((v - 0.1)**2 + 1) + v - 0.1)

        return v_result

    def Ld_difraction(self,d,f,a,deltaN,hrs,hts,hte,hre,hsr,hst,h0,hn,di,hi,omega): # -*- coding: utf-8 -*-

        kbeta = 3
        polarization = 'horizontal'
        Hi = [0,0,0]
        aobt_1 = [-np.inf,-np.inf,-np.inf]
        aobr_1 = [-np.inf,-np.inf,-np.inf]
        Stim_max = -np.inf
        Srim_max = -np.inf
        vmax = -np.inf
        vb_max = -np.inf

        k50 = 157/(157 - deltaN)
        ae = 6371*k50
        abeta = 6371*kbeta

        Ce = 1/ae
        lamb =(3*(10**8))/f



        ##To calculate Lbulls
        val_hi =len(hi)
        for i in range(0, val_hi,1):
            Hi[i] = hi[i] - (hts*(d - di[i]) + hrs*di[i])/d
            aobt_1[i] = (Hi[i]/di[i])

            index = np.where(di[i]-d == 0)
            if (di[i]-d == 0):
                aobr_1[i] = Hi[i]
            else:

                aobr_1[i] = (Hi[i]/(di[i]-d))

        hobs = max(Hi)
        aobt = max(aobt_1)
        aobr = max(aobr_1)

        gt = aobt/(aobt +aobr)
        gr = aobr/(aobt +aobr)

        if(hobs<=0):
            hstp = hst
            hsrp = hsr
        else:
            hstp = hst - hobs*gt
            hsrp = hsr - hobs*gr

        if(hstp>h0):
            hstd = h0
        else:
            hstd = hstp

        if(hsrp>hn):
            hsrd=hn
        else:
            hsrd=hsrp

        hts_n = hts -  hstd
        hrs_n = hrs -  hsrd

        ## The Bullington part of the diffraction calculation
        for p in range(0, 2 ,1):

            if(p==1):
                hts = hts_n
                hrs = hrs_n
                hi = [0,0,0]

            for i in range(0, val_hi,1):

                Stim = (hi[i] + 500*Ce*di[i]*(d - di[i]) - hts)/di[i]

                if Stim>Stim_max:
                   Stim_max = Stim
            Str = (hrs - hts)/d

            if (Stim_max<Str):  #the path is LoS.
               for i in range(0, val_hi,1):

                   if(d-di[i]==0):
                        vmax_i = (hi[i] + 500*Ce*di[i] - ((hts + hrs*di[i])/d))*(np.sqrt((0.002*d)/(lamb*di[i])))
                   else:

                        vmax_i = (hi[i] + 500*Ce*di[i]*(d-di[i])- ((hts*(d - di[i]) + hrs*di[i])/d))*(np.sqrt((0.002*d)/(lamb*di[i]*(d - di[i]))))

                   if vmax_i>vmax:
                       vmax = vmax_i

               Luc = self.Jfunction(vmax)

            if(Stim_max>=Str):
               for i in range(0, val_hi,1):

                   Srim_max = ((hi[i] + 500*Ce*di[i]*(d - di[i]) - hrs)/(d - di[i]))
                   Dbp = (hrs - hts + Srim_max*d)/Stim_max +Srim_max

                   if (d<1):
                      vb = (hts + Stim_max*Dbp - ((hts*(d - Dbp) + hrs*Dbp)/d))*(np.sqrt((0.002)/(lamb*Dbp*(d - Dbp))))
                   else:
                      vb = (hts + Stim_max*Dbp - ((hts*(d - Dbp) + hrs*Dbp)/d))*(np.sqrt((0.002*d)/(lamb*Dbp*(d - Dbp))))

                   if vb>vb_max:
                      vb_max = vb
               Luc = self.Jfunction(vb_max)

            Lbull = Luc + (1 - np.exp(- Luc/6))*(10 + 0.02*d)

            if(p==0):
               Lbulla = Lbull
            if(p==1):
               Lbulls = Lbull


         #Implementing Ldsph
        hte = hts_n
        hre = hrs_n

        if (a=='p'):
           ap = ae
        if (a=='b'):
           ap = abeta

        m = (250*(d**2))/(ap*(hte + hre))
        c = (hte - hre) / (hte + hre)
        b = 2*np.sqrt((m+1)/(3*m))*np.cos(math.pi/3 + (1/3)*np.arccos((3*c/2)*np.sqrt((3*m)/((m+1)**3))))
        dse1 = ((d/2)*(1+b))
        dse2 = (d - dse1)

        dlos = (np.sqrt(2*ap)*(np.sqrt(0.001*hte) + np.sqrt(0.001*hre)))
        hse = ((hte - 500*(dse1**2)/ap)*dse2 + (hre - 500*(dse2**2)/ap)*dse1)/d
        hreq = 17.456*np.sqrt((dse1*dse2*lamb)/d)

        aem = 500*(d/(np.sqrt(hte)+np.sqrt(hre)))**2


        if (d>=dlos):
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

            Kh = 0.036*(adft*f)**(-1/3)*((er - 1)**2 + (18*sig/f)**2)**(-1/4)
            Kv = Kh*(er**2 + (18*sig/f)**2)**0.5

            if (polarization == 'horizontal'):  #horizontal polarization
               K = Kh
            else:
               K = Kv

            Bdft = (1 + 1.6*K**2 + 0.67*K**4)/(1 + 4.5*K**2 + 1.53*K**4)
            Yt = (0.9575*Bdft*((f**2)/adft)**(1/3))*hte
            Yr = (0.9575*Bdft*((f**2)/adft)**(1/3))*hre
            X = (21.88*Bdft*(f/(adft**2))**(1/3))*d

            if (X>=1.6):
               Fx = 11 + 10*np.log10(X) - 17.6*X
            else:
               Fx = -20*np.log10(X) - 5.6488*(X**1.425)

            if i == 0:
               Ldftland = -Fx - self.func_Gt(K,Bdft,Yt) - self.func_Gr(K,Bdft,Yr)
            if i == 1:
               Ldftsea = -Fx - self.func_Gt(K,Bdft,Yt) - self.func_Gr(K,Bdft,Yr)


        Ldft  = omega*Ldftsea + (1- omega)*Ldftland

        if (d>=dlos):
            Ldsph = Ldft

        if(hse>hreq) or (Ldft<0):
            Ldsph = 0
        else:
            Ldsph = (1 - hse/hreq)*Ldft


        Ld = Lbulla + max(Ldsph - Lbulls, 0)


        return Ld

    def get_loss(self, beta, *args, **kwargs) -> np.array:
         d = np.asarray(kwargs["distance"])*(1e-3)  #Km
         f = np.asarray(kwargs["frequency"])*(1e-3)  #GHz
         Ph = np.asarray(kwargs["atmospheric_pressure"])
         T = np.asarray(kwargs["air_temperature"])
         ro = np.asarray(kwargs["water_vapour"])
         deltaN = np.asarray(kwargs["delta_N"])
         hrs = np.asarray(kwargs["Hrs"])
         hts = np.asarray(kwargs["Hts"])
         hte = np.asarray(kwargs["Hte"])
         hre = np.asarray(kwargs["Hre"])
         hsr = np.asarray(kwargs["Hsr"])
         hst = np.asarray(kwargs["Hst"])
         h0 = np.asarray(kwargs["H0"])
         hn = np.asarray(kwargs["Hn"])
         di = np.asarray(kwargs["dist_di"])
         hi = np.asarray(kwargs["hight_hi"])
         omega = np.asarray(kwargs["omega"])
         dlt = np.asarray(kwargs["Dlt"])
         dlr = np.asarray(kwargs["Dlr"])
         p = np.asarray(kwargs["percentage_p"])
         Beta = beta

         C0 = np.asarray(kwargs["C0"])
         C1 = np.asarray(kwargs["C1"])
         C2 = np.asarray(kwargs["C2"])
         D1 = np.asarray(kwargs["D1"])
         D2 = np.asarray(kwargs["D2"])
         D3 = np.asarray(kwargs["D3"])

         if (p>50):
            if (Beta>p):
                #Approximation to the inverse cumulative normal distribution function
                Tx1 = np.sqrt((-2*np.log(p/100)))
                Ex1 = (((C2*Tx1+C1)*Tx1) +C0)/(((D3*Tx1 + D2)*Tx1 + D1)*Tx1 + 1)

                Tx2 = np.sqrt((-2*np.log(beta/100)))
                Ex2 = (((C2*Tx2+C1)*Tx2) +C0)/(((D3*Tx2 + D2)*Tx2 + D1)*Tx2 + 1)

                I1 = Ex1 - Tx1
                I2 = Ex2 - Tx2
                Fi = I1/ I2
            else:
                Fi =1
         else:
            Fi =1

         Ag = self.propagationAg.get_loss(distance=d, frequency=f,atmospheric_pressure=Ph, air_temperature=T, water_vapour=ro)

         Lbfgs = 92.5 + 20*np.log10(f) + 20*np.log10(d) +Ag
         Esp = 2.6*(1 - np.exp(-0.1*(dlr+dlt)))*np.log10(p/50)
         Lb0p = Lbfgs + Esp

         a = 'p'
         Ld50 = self.Ld_difraction(d,f,a,deltaN,hrs,hts,hte,hre,hsr,hst,h0,hn,di,hi,omega)    # Ld50

         a = 'b'
         Ldbeta = self.Ld_difraction(d,f,a,deltaN,hrs,hts,hte,hre,hsr,hst,h0,hn,di,hi,omega)  #Ldbeta

         Ldp = Ld50 + Fi*(Ldbeta - Ld50)
         Ldb = Lb0p + Ldp


         return Ld50, Ldbeta,Ldp, Ldb
         #return Ldb
